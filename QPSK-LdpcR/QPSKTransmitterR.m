% QPSK TransmitterR intro LDPC
classdef QPSKTransmitterR < matlab.System
    properties (Nontunable)
        UpsamplingFactor = 4;
        MessageLength = 105;
        DataLength = 174;
        TransmitterFilterCoefficients = 1;
        ScramblerBase = 2;
        ScramblerPolynomial = [1 1 1 0 1];
        ScramblerInitialConditions = [0 0 0 0];

        % LDPC
        LdpcNewH = zeros(148,296);
        LdpcU = zeros(148,148);
        LdpcL = zeros(148,148);
    end

    properties (Access=private)
        pBitGenerator
        pQPSKModulator 
        pTransmitterFilter        
    end

    methods 
        function obj = QPSKTransmitterR(varargin)
            setProperties(obj, nargin, varargin{:});
        end
    end

    methods (Access=protected)
        function setupImpl(obj)
            obj.pBitGenerator = QPSKBitsGeneratorR( ...
                'MessageLength',                obj.MessageLength, ...
                'BernoulliLength',              obj.DataLength-2*obj.MessageLength,...
                'ScramblerBase',                obj.ScramblerBase, ...
                'ScramblerPolynomial',          obj.ScramblerPolynomial, ...
                'ScramblerInitialConditions',   obj.ScramblerInitialConditions,...
                'LdpcNewH', obj.LdpcNewH,...                % LDPC parameter
                'LdpcU', obj.LdpcU,...
                'LdpcL', obj.LdpcL);
            obj.pQPSKModulator  = comm.QPSKModulator( ...
                'BitInput',                     true, ...
                'PhaseOffset',                  pi/4);
            obj.pTransmitterFilter = dsp.FIRInterpolator( ...
                obj.RolloffFactor, ...
                obj.TransmitterFilterCoefficients);
        end
    end

    methods (Access=protected)
        function [transmittedSignal, transmittedData, modulatedData] = stepImpl(obj)
            % bit to transmitted
            [transmittedData, ~] = step(obj.pBitGenerator);
            
            % modulated bites to QPSK
            modulatedData = step(obj.pQPSKModulator, transmittedData);

            % square root raised cosine 
            transmittedSignal = step(obj.pTransmitterFilter, modulatedData);
        end

        function resetImpl(obj)
            reset(obj.pBitGenerator);
            reset(obj.pQPSKModulator);
            reset(obj.pTransmitterFilter);
        end

        function releaseImpl(obj)
            release(obj.pBitGenerator);
            release(obj.pQPSKModulator);
            release(obj.pTransmitterFilter);
        end

        function N = getNumInputsImpl(~)
            N = 0;
        end 
    end
end
