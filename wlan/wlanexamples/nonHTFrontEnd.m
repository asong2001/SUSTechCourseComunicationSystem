classdef nonHTFrontEnd < matlab.System

% Copyright 2015 The MathWorks, Inc.

%#codegen

properties (Nontunable)
    ChannelBandwidth = 'CBW20';
    SymbolTimingThreshold = .3;
    OFDMSymbolOffset = .75;
end

properties (Access = private, Nontunable)
    pSymbolLen
    pAGC
    pCoarseFreqCompensator
    pFineFreqCompensator
    pSyncSymbolBuffer
end

properties (Access = private)    
    % Packet detection related 
    pPacketDetected
    pLSTFSearchBuffer
    pTimingSynced
    % Symbol timing related
    pLSIGDecoded
    pLLTFSearchBuffer
    pLLTFBufferedSymbols
    pTimingOffset
    % CFO related
    pCoarseCFOEst
    pFineCFOEst
    % L-SIG decoding related
    pCfgRec
    % PSDU decoding related
    pNumDataSymbols
    pFullPayload
    pNumCollectedDataSym
    % Output related
    pChanEst
    pCfgNonHT
    pNoiseVarEst
end

properties (Constant, Hidden)
  ChannelBandwidthSet = matlab.system.StringSet({'CBW20', 'CBW40', 'CBW80', 'CBW160'});
end

methods
  function obj = nonHTFrontEnd(varargin)
    setProperties(obj,nargin,varargin{:});
  end
  
  function set.SymbolTimingThreshold(obj, val)
    prop = 'SymbolTimingThreshold';
    validateattributes(val, {'double'}, ...
        {'real','scalar','>',0,'<',1}, ...
        [class(obj) '.' prop], prop);
    obj.(prop) = val;
  end
  
  function set.OFDMSymbolOffset(obj, val)
    prop = 'OFDMSymbolOffset';
    validateattributes(val, {'double'}, ...
        {'real','scalar','>=',0,'<=',1}, ...
        [class(obj) '.' prop], prop);
    obj.(prop) = val;
  end
  
end

methods (Access = protected)
  function validateInputsImpl(~, x)
    validateattributes(x, {'double'}, {'column','finite'}, 'step', 'X');
  end

  function setupImpl(obj)
    Rs = real(helperSampleRate(obj.ChannelBandwidth));
    obj.pSymbolLen = 4*Rs/1e6;

    % Instantiate objects
    obj.pAGC = comm.AGC;

    obj.pCoarseFreqCompensator = comm.PhaseFrequencyOffset( ...
        'FrequencyOffsetSource', 'Input port', ...
        'SampleRate',            Rs);

    obj.pFineFreqCompensator = comm.PhaseFrequencyOffset( ...
        'FrequencyOffsetSource', 'Input port', ...
        'SampleRate',            Rs);

    obj.pSyncSymbolBuffer = dsp.VariableIntegerDelay( ... 
        'MaximumDelay',      obj.pSymbolLen, ...
        'DirectFeedthrough', true);

    obj.pCfgRec = wlanRecoveryConfig( ...
        'OFDMSymbolOffset',   obj.OFDMSymbolOffset, ...
        'EqualizationMethod', 'ZF'); % Set to zero-forcing
    
    obj.pCfgNonHT = wlanNonHTConfig;
  end

  function resetImpl(obj)
    % Initialize some scalar variables
    obj.pLLTFBufferedSymbols = 0;
    obj.pTimingOffset        = 0;
    obj.pCoarseCFOEst        = 0;
    obj.pFineCFOEst          = 0;
    obj.pNumDataSymbols      = 0;
    obj.pNumCollectedDataSym = 0;
    obj.pNoiseVarEst         = 0; 
    
    % Initialize flags for modes
    obj.pPacketDetected = false;
    obj.pTimingSynced   = false;
    obj.pLSIGDecoded    = false;
      
    % Initialize buffers and states
    obj.pLSTFSearchBuffer = complex(zeros(2*obj.pSymbolLen, 1));
    obj.pLLTFSearchBuffer = complex(zeros(4*obj.pSymbolLen, 1));  
    obj.pFullPayload      = complex(0); 
    obj.pChanEst          = complex(0);

    % Reset System objects
    reset(obj.pAGC);
    reset(obj.pCoarseFreqCompensator);
    reset(obj.pFineFreqCompensator);
    reset(obj.pSyncSymbolBuffer);
  end
  
  function [validPacket, cfgNonHT, rxNonHTData, chanEst, noiseVarEst, data] = stepImpl(obj, x)    
    % Output initialization
    validPacket = false;
    cfgNonHT    = obj.pCfgNonHT;
    rxNonHTData = complex(0);
    chanEst     = complex(0);
    noiseVarEst = 0;    
    
    % Parameters
    chanBW = obj.ChannelBandwidth;
    symLen = obj.pSymbolLen;
    numSym = floor(length(x)/symLen); % Number of symbols in input
    
    for symIdx = 1:numSym % Process input symbol-by-symbol
        data = x((symIdx-1)*symLen + (1:symLen));
        
        % Keep updating L-STF search buffer in case of early failure
        obj.pLSTFSearchBuffer = [obj.pLSTFSearchBuffer(symLen+1:end); data];
        
        if ~obj.pPacketDetected % Packet Detect
            packetStartIdx = helperPacketDetect(obj.pLSTFSearchBuffer, chanBW);
            
            if ~isempty(packetStartIdx) && (packetStartIdx(1) <= symLen)
                % Estimate CFO when more than one L-STF symbol in buffer
                obj.pCoarseCFOEst = real(wlanCoarseCFOEstimate( ...
                    obj.pLSTFSearchBuffer(packetStartIdx(1):end,:), chanBW));
                
                % Switch to symbol timing mode
                obj.pPacketDetected     = true;                
                obj.pTimingSynced       = false;
                obj.pLLTFSearchBuffer   = complex(zeros(size(obj.pLLTFSearchBuffer)));
                obj.pLLTFBufferedSymbols = 0;
            end
        else
            % AGC
            data = step(obj.pAGC, data);
            
            % Coarse frequency offset compensator
            data = step(obj.pCoarseFreqCompensator, data, -obj.pCoarseCFOEst);
            
            if ~obj.pTimingSynced % Symbol timing
                % Update L-LTF search buffer
                obj.pLLTFBufferedSymbols = obj.pLLTFBufferedSymbols + 1;
                obj.pLLTFSearchBuffer((obj.pLLTFBufferedSymbols-1)*symLen ...
                    + (1:symLen), :) = data;
                LLTFStart = helperSymbolTiming(obj.pLLTFSearchBuffer, ...
                    chanBW, obj.SymbolTimingThreshold);

                LLTFLen = 2*symLen;
                % L-LTF Found & the whole L-LTF is in the buffer
                if ~isempty(LLTFStart) && (LLTFStart(1) > 0) && ...
                        (obj.pLLTFBufferedSymbols*symLen - LLTFStart(1) + 1 >= LLTFLen)

                    % Extract L-LTF
                    LLTF = obj.pLLTFSearchBuffer(LLTFStart(1):LLTFStart(1)+LLTFLen-1);

                    % Fine frequency offset compensator
                    obj.pFineCFOEst = real(wlanFineCFOEstimate(LLTF, chanBW));
                    LLTF(1:end/2) = step(obj.pFineFreqCompensator, ...
                        LLTF(1:symLen), -obj.pFineCFOEst);
                    LLTF(end/2+1:end) = step(obj.pFineFreqCompensator, ...
                        LLTF(symLen+(1:symLen)), -obj.pFineCFOEst);

                    % Channel estimation
                    demodLLTF = wlanLLTFDemodulate(LLTF, chanBW);
                    obj.pChanEst = wlanLLTFChannelEstimate(demodLLTF, chanBW);
                    
                    % Estimate noise power using L-LTF field
                    obj.pNoiseVarEst = helperNoiseEstimate(demodLLTF);
                    
                    % Extract L-SIG samples, if any, from L-LTF saerch buffer
                    leftLSIGSamp = obj.pLLTFSearchBuffer(LLTFStart(1)+LLTFLen:obj.pLLTFBufferedSymbols*symLen);
                    obj.pTimingOffset = size(leftLSIGSamp, 1);

                    % Perform symbol synchronization
                    symSyncInput = [complex(zeros(symLen-obj.pTimingOffset, 1)); ...
                        leftLSIGSamp(1:obj.pTimingOffset,:)];
                    step(obj.pSyncSymbolBuffer, symSyncInput(1:symLen,:), obj.pTimingOffset);                        

                    % Switch to L-SIG decoding.
                    obj.pTimingSynced = true;                
                    obj.pLSIGDecoded = false;
                elseif obj.pLLTFBufferedSymbols == 4 
                    % Symbol timing failed -- switch back to packet detection 
                    obj.pPacketDetected = false;
                    obj.pTimingSynced = false;
                end
            else % L-SIG decoding and PSDU buffering
                % Perform symbol synchronization
                syncedSym = step(obj.pSyncSymbolBuffer, complex(data), obj.pTimingOffset);

                % Fine frequency offset compensator
                syncedSym(1:symLen,:) = step(obj.pFineFreqCompensator, syncedSym(1:symLen,:), -obj.pFineCFOEst);

                if ~obj.pLSIGDecoded % L-SIG decoding
                    [LSIGBits, failParityCheck] = wlanLSIGRecover(...
                        syncedSym, obj.pChanEst, obj.pNoiseVarEst, chanBW, obj.pCfgRec);
                    
                    % L-SIG evaluation
                    if ~failParityCheck
                        % Recover packet parameters
                        rate = bi2de(double(LSIGBits(1:3).'),  'left-msb');
                        if rate <= 1
                            obj.pCfgNonHT.MCS = rate + 6;
                        else
                            obj.pCfgNonHT.MCS = mod(rate, 6);
                        end                    
                        obj.pCfgNonHT.PSDULength = bi2de(double(LSIGBits(6:17)'));

                        % Obtain number of OFDM symbols in data field
                        obj.pNumDataSymbols = getNumDataSymbols(obj);
                        
                        % Switch to PSDU buffering mode
                        obj.pLSIGDecoded = true;
                        obj.pFullPayload = complex(zeros(obj.pNumDataSymbols*symLen, 1));
                        obj.pNumCollectedDataSym = 0;
                    else % L-SIG parity failed -- switch back to packet detection 
                        obj.pPacketDetected = false;
                    end
                else % PSDU buffering
                    % Keep buffering payload
                    obj.pNumCollectedDataSym = obj.pNumCollectedDataSym + 1;
                    obj.pFullPayload((obj.pNumCollectedDataSym-1)*symLen+(1:symLen), :) = syncedSym(1:symLen, :);

                    if obj.pNumCollectedDataSym == obj.pNumDataSymbols
                        % Output when payload is full
                        validPacket = true;
                        cfgNonHT    = obj.pCfgNonHT;
                        rxNonHTData = obj.pFullPayload(1:obj.pNumDataSymbols*symLen, :);
                        chanEst     = obj.pChanEst;
                        noiseVarEst = obj.pNoiseVarEst;
                        
                        % Switch back to packet detection
                        obj.pPacketDetected = false;
                    end
                end
            end 
        end 
    end 
  end
  
  function flag = isInputComplexityLockedImpl(~,~)
    flag = false;
  end
  
  function releaseImpl(obj)
    % Release System objects
    release(obj.pAGC);
    release(obj.pCoarseFreqCompensator);
    release(obj.pFineFreqCompensator);
    release(obj.pSyncSymbolBuffer);
  end
  
end

methods (Access = private)
    function numDataSym = getNumDataSymbols(obj)
    % Get number of OFDM data symbols
    
        numSD = 48;       % Data subcarriers
        switch obj.pCfgNonHT.MCS
          case 0 % 6 Mbps
            numBPSCS = 1;  % 'BPSK'
            rate = 1/2;
          case 1 % 9 Mbps
            numBPSCS = 1; 
            rate   = 3/4;
          case 2 % 12 Mbps
            numBPSCS = 2;  % QPSK
            rate   = 1/2;
          case 3 % 18 Mbps
            numBPSCS = 2; 
            rate   = 3/4;
          case 4 % 24 Mbps
            numBPSCS = 4;  % 16QAM 
            rate   = 1/2;
          case 5 % 36 Mbps
            numBPSCS = 4;  
            rate   = 3/4;
          case 6  % 48 Mbps
            numBPSCS = 6;  % '64QAM'
            rate   = 2/3;
          otherwise % 7 => 54 Mbps
            numBPSCS = 6;
            rate   = 3/4;
        end    
        numCBPS = numSD * numBPSCS;
        numDBPS = numCBPS * rate;          
        Ntail = 6; Nservice = 16;
        numDataSym = ceil((8*obj.pCfgNonHT.PSDULength + Nservice + Ntail)/numDBPS);
    end        
end

end

% [EOF]