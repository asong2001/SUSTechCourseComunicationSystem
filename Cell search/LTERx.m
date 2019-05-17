%% Receiver Design: System Architecture
% The general structure of the LTE receiver can be described as follows:
%
% # 1. 首先捕获信号；
% # 2. 信号频偏纠正
% # 3. 信号同步
% # 4. OFDM解调
% # 5. 信道估计
% # 6. 解码PDSCH and DL-SCH 
% # 7. 图像重构

clc
clear
%1. 捕获USRP信号
%(1)设置USRP参数；注意设置要参考发射机参数，特别是帧大小，即发射波形点数。
prmQPSKReceiver.USRPCenterFrequency = 900e6;
prmQPSKReceiver.USRPGain = 25;
prmQPSKReceiver.RxBufferedFrames=1;
prmQPSKReceiver.Fs = 5e6; % IQ 速率；
prmQPSKReceiver.USRPDecimationFactor = 100e6/prmQPSKReceiver.Fs; % USRP上采样因子是100M/20M=5
prmQPSKReceiver.FrameSize=307200;
prmQPSKReceiver.USRPFrameLength = prmQPSKReceiver.FrameSize*prmQPSKReceiver.RxBufferedFrames;

%(2) 构造接收机对象
radio = comm.SDRuReceiver(...
        'IPAddress',            '192.168.10.2', ...
        'CenterFrequency',      prmQPSKReceiver.USRPCenterFrequency, ...
        'Gain',                 prmQPSKReceiver.USRPGain, ...
        'DecimationFactor',     prmQPSKReceiver.USRPDecimationFactor, ...
        'FrameLength',          prmQPSKReceiver.USRPFrameLength, ...
        'OutputDataType',       'double');

%(3) 循环捕获，直到成功捕获数据包    
errorIndex=0;
while (true)
   % 从USRP读取IQ信号；
    [corruptSignal, len] = step(radio);

   % 如果未能成功读取我们希望的长度，报错
      if len < prmQPSKReceiver.USRPFrameLength 
         errorIndex = errorIndex+1;
         disp ( 'Not enough samples returned!' ) ;
         disp(errorIndex)
      else
          
         rxWaveform1 = corruptSignal;    
         break; % 否则，跳出循环，进一步做数据包解码操作。 
      end
end

%%
% *Receiver Setup*

% User defined parameters --- configure the same as transmitter
rxsim = struct;
rxsim.RadioFrontEndSampleRate = 15.36e6; % Configure for same sample rate
                                                       % as transmitter
rxsim.RadioCenterFrequency = 900e6; % Center frequency in Hz
rxsim.NRxAnts = 1;      % Number of transmit antennas
rxsim.FramesPerBurst = 2+1; % Number of LTE frames to capture in each burst.
                                          % Capture 1 more LTE frame than transmitted to  
                                          % allow for timing offset wraparound...
rxsim.numBurstCaptures = 1; % Number of bursts to capture

% Derived parameters
samplesPerFrame = 10e-3*rxsim.RadioFrontEndSampleRate; % LTE frames period is 10 ms
% samplesPerFrame = 768000;

%%

rx.BasebandSampleRate = rxsim.RadioFrontEndSampleRate;
rx.CenterFrequency = rxsim.RadioCenterFrequency;
rx.SamplesPerFrame = samplesPerFrame;
rx.OutputDataType = 'double';
rx.EnableBurstMode = true;
rx.NumFramesInBurst = rxsim.FramesPerBurst;

rx.ChannelMapping = 1;

% burstCaptures holds rx.FramesPerBurst number of consecutive frames worth
% of baseband LTE samples. Each column holds one LTE frame worth of data.
burstCaptures = zeros(samplesPerFrame,rxsim.NRxAnts,rxsim.FramesPerBurst);

%%
% *LTE Receiver Setup* 
txsim.RC = 'R.7';
rmc = lteRMCDL(txsim.RC);

enb.PDSCH = rmc.PDSCH;
enb.DuplexMode = 'FDD';
enb.CyclicPrefix = 'Normal';
enb.CellRefP = 4; 

%%
% Bandwidth: {1.4 MHz, 3 MHz, 5 MHz, 10 MHz, 20 MHz}
SampleRateLUT = [1.92 3.84 7.68 15.36 30.72]*1e6;
NDLRBLUT = [6 15 25 50 100];
enb.NDLRB = NDLRBLUT(SampleRateLUT==rxsim.RadioFrontEndSampleRate);
if isempty(enb.NDLRB)
    error('Sampling rate not supported. Supported rates are %s.',...
            '1.92 MHz, 3.84 MHz, 7.68 MHz, 15.36 MHz, 30.72 MHz');
end
fprintf('\nSDR hardware sampling rate configured to capture %d LTE RBs.\n',enb.NDLRB);

%%
% Channel estimation configuration structure
cec.PilotAverage = 'UserDefined';  % Type of pilot symbol averaging
cec.FreqWindow = 9;                % Frequency window size in REs
cec.TimeWindow = 9;                % Time window size in REs
cec.InterpType = 'Cubic';          % 2D interpolation type
cec.InterpWindow = 'Centered';     % Interpolation window type
cec.InterpWinSize = 3;             % Interpolation window size

enbDefault = enb;
%%
% *Signal Capture and Processing*

%    rxWaveform  = eNodeBOutput;
    
    % Show power spectral density of captured burst
%    hsa.SampleRate = rxsim.RadioFrontEndSampleRate;
%    step(hsa,rxWaveform);

%      load('eNodeBOutput.mat');    
%      rxWaveform  = eNodeBOutput;
%      
%    load('rxWaveform1.mat');
    rxWaveform  = [rxWaveform1;rxWaveform1];

    % Perform frequency offset correction for known cell ID
    frequencyOffset = lteFrequencyOffset(enb,rxWaveform);
    rxWaveform = lteFrequencyCorrect(enb,rxWaveform,frequencyOffset);
    fprintf('\nCorrected a frequency offset of %i Hz.\n',frequencyOffset)
    
    % Perform the blind cell search to obtain cell identity and timing offset
    %   Use 'PostFFT' SSS detection method to improve speed
    cellSearch.SSSDetection = 'PostFFT'; cellSearch.MaxCellCount = 1;
    [NCellID,frameOffset] = lteCellSearch(enb,rxWaveform,cellSearch);
    fprintf('Detected a cell identity of %i.\n', NCellID);
    enb.NCellID = NCellID; % From lteCellSearch
    
    % Sync the captured samples to the start of an LTE frame, and trim off
    % any samples that are part of an incomplete frame.
    rxWaveform = rxWaveform(frameOffset+1:end,:);
    tailSamples = mod(length(rxWaveform),samplesPerFrame);
    rxWaveform = rxWaveform(1:end-tailSamples,:);
    enb.NSubframe = 0;
    fprintf('Corrected a timing offset of %i samples.\n',frameOffset)
    
    % OFDM demodulation
    rxGrid = lteOFDMDemodulate(enb,rxWaveform);
    
    % Perform channel estimation for 4 CellRefP as currently we do not
    % know the CellRefP for the eNodeB.
    [hest,nest] = lteDLChannelEstimate(enb,cec,rxGrid);
    
    sfDims = lteResourceGridSize(enb);
    Lsf = sfDims(2); % OFDM symbols per subframe
    LFrame = 10*Lsf; % OFDM symbols per frame
    numFullFrames = length(rxWaveform)/samplesPerFrame;
    
    rxDataFrame = zeros(sum(enb.PDSCH.TrBlkSizes(:)),numFullFrames);
    recFrames = zeros(numFullFrames,1);
    rxSymbols = []; txSymbols = [];
    
    % For each frame decode the MIB, PDSCH and DL-SCH
    for frame = 0:(numFullFrames-1)
        fprintf('\nPerforming DL-SCH Decode for frame %i of %i in burst:\n', ...
            frame+1,numFullFrames)
        
        % Extract subframe #0 from each frame of the received resource grid
        % and channel estimate.
        enb.NSubframe = 0;
        rxsf = rxGrid(:,frame*LFrame+(1:Lsf),:);
        hestsf = hest(:,frame*LFrame+(1:Lsf),:,:);
               
        % PBCH demodulation. Extract resource elements (REs)
        % corresponding to the PBCH from the received grid and channel
        % estimate grid for demodulation.
        enb.CellRefP = 4;
        pbchIndices = ltePBCHIndices(enb); 
        [pbchRx,pbchHest] = lteExtractResources(pbchIndices,rxsf,hestsf);
        [~,~,nfmod4,mib,CellRefP] = ltePBCHDecode(enb,pbchRx,pbchHest,nest);
        
        % If PBCH decoding successful CellRefP~=0 then update info
        if ~CellRefP
            fprintf('  No PBCH detected for frame.\n');
            continue;
        end
        enb.CellRefP = CellRefP; % From ltePBCHDecode
        
        % Decode the MIB to get current frame number
        enb = lteMIB(mib,enb);

        % Incorporate the nfmod4 value output from the function
        % ltePBCHDecode, as the NFrame value established from the MIB
        % is the system frame number modulo 4.
        enb.NFrame = enb.NFrame+nfmod4;
        fprintf('  Successful MIB Decode.\n')
        fprintf('  Frame number: %d.\n',enb.NFrame);
        
        % The eNodeB transmission bandwidth may be greater than the
        % captured bandwidth, so limit the bandwidth for processing
        enb.NDLRB = min(enbDefault.NDLRB,enb.NDLRB);
        
        % Store received frame number
        recFrames(frame+1) = enb.NFrame;
               
        % Process subframes within frame (ignoring subframe 5)
        for sf = 0:9
            if sf~=5 % Ignore subframe 5
                % Extract subframe
                enb.NSubframe = sf;
                rxsf = rxGrid(:,frame*LFrame+sf*Lsf+(1:Lsf),:);

                % Perform channel estimation with the correct number of CellRefP
                [hestsf,nestsf] = lteDLChannelEstimate(enb,cec,rxsf);

                % PCFICH demodulation. Extract REs corresponding to the PCFICH
                % from the received grid and channel estimate for demodulation.
                pcfichIndices = ltePCFICHIndices(enb);
                [pcfichRx,pcfichHest] = lteExtractResources(pcfichIndices,rxsf,hestsf);
                [cfiBits,recsym] = ltePCFICHDecode(enb,pcfichRx,pcfichHest,nestsf);

                % CFI decoding
                enb.CFI = lteCFIDecode(cfiBits);
                
                % Get PDSCH indices
                [pdschIndices,pdschIndicesInfo] = ltePDSCHIndices(enb, enb.PDSCH, enb.PDSCH.PRBSet); 
                [pdschRx, pdschHest] = lteExtractResources(pdschIndices, rxsf, hestsf);

                % Perform deprecoding, layer demapping, demodulation and
                % descrambling on the received data using the estimate of
                % the channel
                [rxEncodedBits, rxEncodedSymb] = ltePDSCHDecode(enb,enb.PDSCH,pdschRx,...
                                               pdschHest,nestsf);

                % Append decoded symbol to stream
                rxSymbols = [rxSymbols; rxEncodedSymb{:}]; %#ok<AGROW>

                % Transport block sizes
                outLen = enb.PDSCH.TrBlkSizes(enb.NSubframe+1);  

                % Decode DownLink Shared Channel (DL-SCH)
                [decbits{sf+1}, blkcrc(sf+1)] = lteDLSCHDecode(enb,enb.PDSCH,...
                                                outLen, rxEncodedBits);  %#ok<SAGROW>

                % Recode transmitted PDSCH symbols for EVM calculation                            
                %   Encode transmitted DLSCH 
                txRecode = lteDLSCH(enb,enb.PDSCH,pdschIndicesInfo.G,decbits{sf+1});
                %   Modulate transmitted PDSCH
                txRemod = ltePDSCH(enb, enb.PDSCH, txRecode);
                %   Decode transmitted PDSCH
                [~,refSymbols] = ltePDSCHDecode(enb, enb.PDSCH, txRemod);
                %   Add encoded symbol to stream
                txSymbols = [txSymbols; refSymbols{:}]; %#ok<AGROW>

%               release(hcd); % Release previous constellation plot
%               step(hcd,rxEncodedSymb{:}); % Plot current constellation
            end
        end
        
        % Reassemble decoded bits
        fprintf('  Retrieving decoded transport block data.\n');
        rxdata = [];
        for i = 1:length(decbits)
            if i~=6 % Ignore subframe 5
                rxdata = [rxdata; decbits{i}{:}]; %#ok<AGROW>
            end
        end
        
        % Store data from receive frame
        rxDataFrame(:,frame+1) = rxdata;

        % Plot channel estimate between CellRefP 0 and the receive antennae
        focalFrameIdx = frame*LFrame+(1:LFrame);
              
    end

% Release both transmit and receive objects once reception is complete

%%
% *Result Qualification and Display*

% Determine index of first transmitted frame (lowest received frame number)
[~,frameIdx] = min(recFrames);

fprintf('\nRecombining received data blocks:\n');

decodedRxDataStream = zeros(length(rxDataFrame(:)),1);
frameLen = size(rxDataFrame,1);
% Recombine received data blocks (in correct order) into continuous stream
for n=1:numFullFrames
    currFrame = mod(frameIdx-1,numFullFrames)+1; % Get current frame index 
    decodedRxDataStream((n-1)*frameLen+1:n*frameLen) = rxDataFrame(:,currFrame);
    frameIdx = frameIdx+1; % Increment frame index
end

% Perform EVM calculation
if ~isempty(rxSymbols)
    hEVM = comm.EVM();
    hEVM.MaximumEVMOutputPort = true;
    [evm.RMS,evm.Peak] = step(hEVM,txSymbols, rxSymbols);
    fprintf('  EVM peak = %0.3f%%\n',evm.Peak);
    fprintf('  EVM RMS  = %0.3f%%\n',evm.RMS);
else
    fprintf('  No transport blocks decoded.\n');
end

% Perform bit error rate (BER) calculation
% hBER = comm.ErrorRate;
% err = step(hBER, decodedRxDataStream(1:length(trData)), trData);
% fprintf('  Bit Error Rate (BER) = %0.5f.\n', err(1));
% fprintf('  Number of bit errors = %d.\n', err(2));
% fprintf('  Number of transmitted bits = %d.\n',length(trData));

% Recreate image from received data size(decodedRxDataStream)=545888
fprintf('\nConstructing image from received data.\n');
% 只收一个帧可以这样处理
% decodedRxDataStream=[decodedRxDataStream;decodedRxDataStream];
decodedRxDataStream=decodedRxDataStream(1:545888);
str = reshape(sprintf('%d',decodedRxDataStream(1:422280)), 8, []).';
decdata = uint8(bin2dec(str));
imsize =[115   153     3];
receivedImage = reshape(decdata,imsize);

% Plot receive image

figure(1);
subplot(212); 
imshow(receivedImage);
%title(sprintf('Received Image: %dx%d Antenna Configuration',txsim.NTxAnts, rxsim.NRxAnts));
