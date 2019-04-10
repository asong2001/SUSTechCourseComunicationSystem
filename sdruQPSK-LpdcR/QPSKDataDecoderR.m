classdef QPSKDataDecoderR < matlab.System
    
    % Copyright 2012-2017 The MathWorks, Inc.
    
    properties (Nontunable)
        FrameSize=100;
        BarkerLength=13;
        ModulationOrder = 4; 
        DataLength=174;
        MessageLength=112;
%         HeaderLength = 26;
%         PayloadLength = 2240;
%         NumberOfMessage = 20;
        DescramblerBase = 2;
        DescramblerPolynomial = [1 1 1 0 1];
        DescramblerInitialConditions = [0 0 0 0];
%         BerMask = [];
        % ldpc parameter
        LdpcNewH = zeros(148,296);
        LdpcU = zeros(148,148);
        LdpcL = zeros(148,148);
        LdpcIteration = 1;
        PrintOption = false;
    end
    
    properties (Access = private)
        pCorrelator
%         pPayloadLength
        pQPSKDemodulator
        pDescrambler
        pBitGenerator
        pErrorRateCalc
%         pTargetBits
%         pBitToInteger
%         pIntegerToBit
        pBER
    end
    
    properties (Constant, Access = private)
        pBarkerCode = [+1; +1; +1; +1; +1; -1; -1; +1; +1; -1; +1; -1; +1;...
                       +1; +1; +1; +1; +1; -1; -1; +1; +1; -1; +1; -1; +1]; % Bipolar Barker Code
        pModulatedHeader = sqrt(2)/2 * (-1-1i) * QPSKDataDecoderR.pBarkerCode;
%         pMessage = 'Hello world';
%         pMessageLength = 16;
    end
    
    methods
        function obj = QPSKDataDecoderR(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods (Access = protected)
        function setupImpl(obj, ~)
%             coder.extrinsic('sprintf');
            obj.pCorrelator = dsp.Crosscorrelator;
            obj.pQPSKDemodulator = comm.QPSKDemodulator('PhaseOffset',pi/4, ...
                'BitOutput', true);
            
            obj.pDescrambler = comm.Descrambler(obj.DescramblerBase, ...
                obj.DescramblerPolynomial, obj.DescramblerInitialConditions);
            
            obj.pBitGenerator = QPSKBitsGeneratorR('MessageLength',obj.MessageLength,...
                'BernoulliLength',obj.DataLength-2*obj.MessageLength,...
                'ScramblerBase',obj.DescramblerBase,...
                'ScramblerPolynomial',obj.DescramblerPolynomial,...
                'ScramblerInitialConditions',obj.DescramblerInitialConditions,...
                'LdpcNewH',obj.LdpcNewH,...
                'LdpcU',obj.LdpcU,...
                'LdpcL',obj.LdpcL);
            
            
            obj.pErrorRateCalc = comm.ErrorRate;
%             ( ...
%                 'Samples', 'Custom', ...
%                 'CustomSamples', obj.BerMask);
%             
%             obj.pBitToInteger = comm.BitToInteger(7, 'OutputDataType', 'int8');  % 7 bits per character
%             obj.pIntegerToBit = comm.IntegerToBit(7, 'OutputDataType', 'double');% Double for BER calculation
%             
%             % Since we only calculate BER on message part, 000s are not
%             % necessary here, they are just place-holder.
%             msgSet = zeros(obj.NumberOfMessage * obj.pMessageLength, 1);
%             for msgCnt = 0 : obj.NumberOfMessage - 1
%                 msgSet(msgCnt * obj.pMessageLength + (1 : obj.pMessageLength)) = ...
%                     sprintf('%s %03d\n', obj.pMessage, mod(msgCnt, 100));
%             end
%             obj.pTargetBits = obj.pIntegerToBit(msgSet);
        end
        
        function  BER = stepImpl(obj, data)
%             if isValid
                % Phase ambiguity estimation
                phaseEst = round(angle(mean(conj(obj.pModulatedHeader) .* data(1:obj.BarkerLength)))*2/pi)/2*pi;
                
                % Compensating for the phase ambiguity
                phShiftedData = data .* exp(-1i*phaseEst);
                
                % Demodulating the phase recovered data
                demodOut = step(obj.pQPSKDemodulator,phShiftedData);

                demodOutMsg = demodOut(...
                    obj.BarkerLength*log2(obj.ModulationOrder)+1:...
                    obj.FrameSize*log2(obj.ModulationOrder));

                vhat = decodeBitFlip(demodOutMsg', obj.LdpcNewH, obj.LdpcIteration);
                deScrDataMsg = vhat(149:end)';
                
                % Performs descrambling on only payload part
                deScrData = step(obj.pDescrambler,deScrDataMsg);
%                     demodOut(obj.HeaderLength + (1 : obj.PayloadLength)));

                % Recovering the message from the data
                Received = deScrData(1:obj.MessageLength);
                bits2ASCII(obj,Received);
                
                [~,transmittedMessage]=step(obj.pBitGenerator);
                
                BER=step(obj.pErrorRateCalc,transmittedMessage,Received);
        end
%                 if(obj.PrintOption)
%                     fprintf('%s', char(charSet));
%                 end
%                 
%                 % Perform BER calculation only on message part
%                 obj.pBER = obj.pErrorRateCalc(obj.pTargetBits, deScrData);
%             end
%             BER = obj.pBER;
%         end
        
        function resetImpl(obj)
            reset(obj.pCorrelator);
            reset(obj.pQPSKDemodulator);
            reset(obj.pDescrambler);
            reset(obj.pBitGenerator);
            reset(obj.pErrorRateCalc);
%             reset(obj.pIntegerToBit);
%             obj.pBER = zeros(3, 1);
        end
        
        function releaseImpl(obj)
            release(obj.pCorrelator);
            release(obj.pQPSKDemodulator);
            release(obj.pDescrambler);
            release(obj.pBitGenerator);
            release(obj.pErrorRateCalc);
%             release(obj.pIntegerToBit);
%             obj.pBER = zeros(3, 1);
        end
    end
    
    methods(Access=private)
        function bits2ASCII(obj,u)
            coder.extrinsic('disp')
            
            %convert binary-valumed column vector to 7-bit decimal values
            w=[64 32 16 8 4 2 1];
            Nbits=numel(u);
            Ny=Nbits/7;
            y=zeros(1,Ny);
            for i=0:Ny-1
                y(i+1)=w*u(7*i+(1:7));
            end
            
            %display ASCII message to command window
            if(obj.PrintOption)
                disp(char(y));
            end
        end
    end
end

