classdef QPSKReceiverR < matlab.System
    
    % Copyright 2012-2017 The MathWorks, Inc.
    
    properties (Nontunable)
        DesireAmplitude = 1/sqrt(2);
        ModulationOrder = 4;
        DownsamplingFactor = 2;
        CoarseFrequencyResolution = 50;
        PhaseRecoveryLoopBandwidth = 0.01;
        PhaseRecoveryDampingFactor = 1;
        TimingRecoveryDampingFactor = 1;
        TimingRecoveryLoopBandwidth = 0.01;
        TimingErrorDetectorGain = 5.4;
        PostFilterOversampling = 2;
        FrameSize = 100;
        BarkerLength = 26;
        MessageLength = 105;
        SampleRate = 200000;
        DataLength = 148;
        ReceiverCoefficients = 1;
        DescramblerBase = 2;
        DescramblerPolynomial = [1 1 1 0 1];
        DescramblerInitialConditions = [0 0 0 0];
        % ldpc parameter
        LdpcNewH = zeros(148,296);
        LdpcU = zeros(148,148);
        LdpcL = zeros(148,148);
        LdpcIteration = 1;
        PrintOption = false;
    end
    
    properties (Access = private)
        pAGC
        pRxFilter
        pCoarseFreqEstimator
        pCoarseFreqCompensator
        pFineFreqCompensator
        pTimingRec
        pPrbDet
        pFrameSync
        pDataDecod
        pBER
    end
    
    properties (Access = private, Constant)
        pUpdatePeriod = 4 % Defines the size of vector that will be processed in AGC system object
        pBarkerCode = [+1; +1; +1; +1; +1; -1; -1; +1; +1; -1; +1; -1; +1 +1; +1; +1; +1; +1; -1; -1; +1; +1; -1; +1; -1; +1];
        pModulatedHeader = sqrt(2)/2*(-1-1j)*QPSKReceiverR.pBarkerCode;
    end
    
    methods
        function obj = QPSKReceiverR(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods (Access = protected)
        function setupImpl(obj, ~)
            
            obj.pAGC = comm.AGC;
            
            obj.pRxFilter = dps.FIRDecimator( ...
                'Numerator', obj.ReceiverCoefficients,...
                'DecimationFator', obj.DownsamplingFactor);
            
            obj.pCoarseFreqEstimator = comm.CoarseFrequencyCompensator( ...
                'Modulation',               obj.ModulationOrder, ...
                'Algorithm',                'FFT-based', ...
                'MaximumFrequencyOffset',   obj.CoarseFrequencyResolution, ....
                'SampleRate',               obj.SampleRate);
            
            obj.pCoarseFreqCompensator = comm.PhaseFrequencyOffset( ...
                'PhaseOffset',              0, ...
                'FrequencyOffsetSource',    'Input port', ...
                'SampleRate',               obj.SampleRate);

            obj.pMeanFreqOff = 0;
            
            obj.pCnt = 0;
            
            obj.pFineFreqCompensator = comm.CarrierSynchronizer( ...
                'Modulation',               'QPSK', ...
                'ModulationPhaseOffset',    'Auto', ...
                'SamplesPerSymbol',         obj.PostFilterOversampling, ...
                'DampingFactor',            obj.PhaseRecoveryDampingFactor, ...
                'NormalizedLoopBandwidth',  obj.PhaseRecoveryLoopBandwidth);
            
            obj.pTimingRec = comm.SymbolSynchronizer( ...
                'TimingErrorDetector',      'Zero-Crossing (decision-directed)', ...
                'SamplesPerSymbol',         obj.PostFilterOversampling, ...
                'DampingFactor',            obj.TimingRecoveryDampingFactor, ...
                'NormalizedLoopBandwidth',  obj.TimingRecoveryLoopBandwidth, ...
                'DetectorGain',             obj.TimingErrorDetectorGain);
            
            % obj.pPrbDet = comm.PreambleDetector(obj.pModulatedHeader, ...
            %     'Input',                    'Symbol', ...
            %     'Threshold',                obj.PreambleDetectorThreshold);
             
            %  obj.pFrameSync = FrameSynchronizer( ...
            %      'OutputFrameLength',        obj.FrameSize, ...
            %      'PreambleLength',           obj.HeaderLength / 2);
            
            obj.pFrameSync = FrameFormation( ...
                'OutputFrameLength',        obj.FrameSize, ...
                'PerformSynchronization', true,...
                'FrameHeader',           obj.pModulatedHeader);
            
            obj.pDataDecod = QPSKDataDecoder( ...
                'ModulationOrder',          obj.ModulationOrder, ...
                'HeaderLength',             obj.HeaderLength, ...
                'NumberOfMessage',          obj.NumberOfMessage, ...
                'PayloadLength',            obj.PayloadLength, ...
                'DescramblerBase',          obj.DescramblerBase, ...
                'DescramblerPolynomial',    obj.DescramblerPolynomial, ...
                'DescramblerInitialConditions', obj.DescramblerInitialConditions, ...
                'BerMask',                  obj.BerMask, ...
                'PrintOption',              obj.PrintOption,...
                'LdpcNewH', obj.LdpcNewH,...     % LDPC parameter
                'LdpcU', obj.LdpcU,...
                'LdpcL', obj.LdpcL,...
                'LdpcIteration', obj.LdpcIteration);
        end
        
        function [RCRxSignal, timingRecSignal, fineCompSignal, BER] = stepImpl(obj, bufferSignal)
            
            AGCSignal = obj.DesireAmplitude*step(obj.pAGC,bufferSignal);                          % AGC control
            RCRxSignal = step(obj.pRxFilter,AGCSignal);                       % Pass the signal through
                                                                         % Square-Root Raised Cosine Received Filter
            freqOffsetEst = step(obj.pCoarseFreqEstimator, RCRxSignal);   % Coarse frequency offset estimation
            % average coarse frequency offset estimate, so that carrier
            % sync is able to lock/converge

            % freqOffsetEst = (freqOffsetEst + obj.pCnt * obj.pMeanFreqOff)/(obj.pCnt+1);
            % obj.pCnt = obj.pCnt + 1;            % update state
            % obj.pMeanFreqOff = freqOffsetEst;
            
            coarseCompSignal = step(obj.pCoarseFreqCompensator,RCRxSignal,...
                -freqOffsetEst);                                         % Coarse frequency compensation
            % timingRecSignal = obj.pTimingRec(coarseCompSignal);          % Symbol timing recovery
            
            fineCompSignal = step(obj.pFineFreqCompensator,coarseCompSignal);  % Fine frequency compensation
            
            % [prbIdx, dtMt] = obj.pPrbDet(fineCompSignal);                % Detect frame header
            [timingRecSignal,timingRecBuffer] = step(obj.pTimingRec,fineCompSignal);

%             [symFrame, isFrameValid] = obj.pFrameSync(fineCompSignal, ...
%                 prbIdx, dtMt);                                           % Frame synchronization
            [symFrame, isFrameValid] = step(obj.pFrameSync,timingRecSignal);                                           % Frame synchronization
            
            if isFrameValid
                obj.pBER = step(obj.pDataDecod,symFrame);
            end

            BER = obj.pBER;
            
        end
        
        function resetImpl(obj)
            obj.pBER = zeros(3,1);
            reset(obj.pAGC);
            reset(obj.pRxFilter);
            reset(obj.pCoarseFreqEstimator);
            reset(obj.pCoarseFreqCompensator);
            reset(obj.pFineFreqCompensator);
            reset(obj.pTimingRec);
            reset(obj.pTimingRec);
            reset(obj.pFrameSync);
            reset(obj.pDataDecod);
            obj.pMeanFreqOff = 0;
            obj.pCnt = 0;
        end
        
        function releaseImpl(obj)
            release(obj.pAGC);
            release(obj.pRxFilter);
            release(obj.pCoarseFreqEstimator);
            release(obj.pCoarseFreqCompensator);
            release(obj.pFineFreqCompensator);
            release(obj.pTimingRec);
            release(obj.pPrbDet);
            release(obj.pFrameSync);
            release(obj.pDataDecod);
        end
        
        function N = getNumOutputsImpl(~)
            N = 4;
        end
    end
end

