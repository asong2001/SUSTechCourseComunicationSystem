classdef QAMTransmitter < matlab.System  
%#codegen
% Generates the QPSK signal to be transmitted
    
%   Copyright 2012-2016 The MathWorks, Inc.
    
    properties (Nontunable)
        UpsamplingFactor = 4;
        MessageLength = 105;
        DataLength = 174;
        TransmitterFilterCoefficients = 1;
        ScramblerBase = 2;
        ScramblerPolynomial = [1 1 1 0 1];
        ScramblerInitialConditions = [0 0 0 0];
    end
    
     properties (Access=private)
        pBitGenerator
        pQPSKModulator 
        pTransmitterFilter
    end
    
    methods
        function obj = QPSKTransmitter(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods (Access=protected)
        function setupImpl(obj)
            obj.pBitGenerator = QPSKBitsGenerator(...
                'MessageLength', obj.MessageLength, ...
                'BernoulliLength', obj.DataLength-obj.MessageLength, ...
                'ScramblerBase', obj.ScramblerBase, ...
                'ScramblerPolynomial', obj.ScramblerPolynomial, ...
                'ScramblerInitialConditions', obj.ScramblerInitialConditions);
            obj.p16QAMModulator  = comm.RectangularQAMModulator(16,'BitInput',true, ...
                'PhaseOffset', 0, 'NormalizationMethod', ...
                'Average power',...
                 'SymbolMapping', 'Custom', ...
                 'CustomSymbolMapping', [11 10 14 15 9 8 12 13 1 0 4 5 3 2 6 7]);
            obj.pTransmitterFilter = dsp.FIRInterpolator(obj.UpsamplingFactor, ...
                obj.TransmitterFilterCoefficients);
        end
        
        function transmittedSignal = stepImpl(obj)
           
            [transmittedData, ~] = obj.pBitGenerator();                % Generates the data to be transmitted           
            modulatedData = obj.p16QAMModulator(transmittedData);       % Modulates the bits into QPSK symbols           
            transmittedSignal = obj.pTransmitterFilter(modulatedData); % Square root Raised Cosine Transmit Filter
        end
        
        function resetImpl(obj)
            reset(obj.pBitGenerator);
            reset(obj.p16QAMModulator );
            reset(obj.pTransmitterFilter);
        end
        
        function releaseImpl(obj)
            release(obj.pBitGenerator);
            release(obj.p16QAMModulator );
            release(obj.pTransmitterFilter);
        end
        
        function N = getNumInputsImpl(~)
            N = 0;
        end
    end
end

