classdef QPSKDataDecoder < matlab.System

% Copyright 2012-2016 The MathWorks, Inc.
     
    properties (Nontunable)
        FrameSize = 100;
        BarkerLength = 13;
        ModulationOrder = 4;
        DataLength = 174;
        MessageLength = 105;
        DescramblerBase = 2;
        DescramblerPolynomial = [1 1 1 0 1];
        DescramblerInitialConditions = [0 0 0 0];
        PrintOption = false;
        numTx = 2;
        numRx = 2;
    end
    
    properties (Access = private)
        pCorrelator
        pQPSKDemodulator
        pDescrambler
        pBitGenerator
        pErrorRateCalc
        phOSTBCComb
    end
    
    properties (Constant, Access = private)
        pBarkerCode = [+1; +1; +1; +1; +1; -1; -1; +1; +1; -1; +1; -1; +1]; % Bipolar Barker Code        
        pModulatedHeader = sqrt(2)/2 * (-1-1i) * QPSKDataDecoder.pBarkerCode;
    end
    
    methods
        function obj = QPSKDataDecoder(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods (Access = protected)
        function setupImpl(obj, ~)            
            obj.pCorrelator = dsp.Crosscorrelator;
            
            obj.pQPSKDemodulator = comm.QPSKDemodulator('PhaseOffset',pi/4, ...
                'BitOutput', true);
            
            obj.pDescrambler = comm.Descrambler(obj.DescramblerBase, ...
                obj.DescramblerPolynomial, obj.DescramblerInitialConditions);
            
            obj.pBitGenerator = QPSKBitsGenerator('MessageLength', obj.MessageLength, ...
                'BernoulliLength', obj.DataLength-obj.MessageLength, ...
                'ScramblerBase', obj.DescramblerBase, ...
                'ScramblerPolynomial', obj.DescramblerPolynomial, ...
                'ScramblerInitialConditions', obj.DescramblerInitialConditions);
            
            obj.pErrorRateCalc = comm.ErrorRate; 
            
%             obj.phOSTBCComb = comm.OSTBCCombiner(...
%                 'NumTransmitAntennas',obj.numTx,...
%                 'NumReceiveAntennas',obj.numRx);
        end
        
        function  BER = stepImpl(obj, data)            
            % Phase offset estimation
            phaseEst = round(angle(mean(conj(obj.pModulatedHeader) .* data(1:obj.BarkerLength)))*2/pi)/2*pi;
            
            % Compensating for the phase offset
            phShiftedData = data .* exp(-1i*phaseEst);
            
%             combinedData = obj.phOSTBCComb(phShiftedData);
            % Demodulating the phase recovered data
            demodOut = obj.pQPSKDemodulator(phShiftedData);
            
            % Performs descrambling
            deScrData = obj.pDescrambler( ...
                demodOut( ...
                obj.BarkerLength*log2(obj.ModulationOrder)+1 : ...
                obj.FrameSize*log2(obj.ModulationOrder)));
            
            % Recovering the message from the data
            Received = deScrData(1:obj.MessageLength);
            bits2ASCII(obj, Received);
            
            [~, transmittedMessage] = obj.pBitGenerator();
            
            BER = obj.pErrorRateCalc(transmittedMessage, Received);
        end
 
        function resetImpl(obj)
            reset(obj.pCorrelator);
            reset(obj.pQPSKDemodulator);
            reset(obj.pDescrambler);
            reset(obj.pBitGenerator);
            reset(obj.pErrorRateCalc);
        end
        
        function releaseImpl(obj)
            release(obj.pCorrelator);
            release(obj.pQPSKDemodulator);
            release(obj.pDescrambler);
            release(obj.pBitGenerator);
            release(obj.pErrorRateCalc);
        end        
    end
    
    methods (Access=private)
        function bits2ASCII(obj,u)
            coder.extrinsic('disp')
            
            % Convert binary-valued column vector to 7-bit decimal values.
            w = [64 32 16 8 4 2 1]; % binary digit weighting
            Nbits = numel(u);
            Ny = Nbits/7;
            y = zeros(1,Ny);
            for i = 0:Ny-1
                y(i+1) = w*u(7*i+(1:7));
            end
            
            % Display ASCII message to command window   
            if(obj.PrintOption)
                disp(char(y));
            end
        end
    end
end

