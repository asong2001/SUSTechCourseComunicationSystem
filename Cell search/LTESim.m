clc
clear

%%
% Setup Spectrum viewer
hsa = dsp.SpectrumAnalyzer( ...
    'SpectrumType',    'Power density', ...
    'SpectralAverages', 10, ...
    'YLimits',         [-150 -60], ...
    'Title',           'Received Baseband LTE Signal Spectrum', ...
    'YLabel',          'Power spectral density');

% Setup the constellation diagram viewer for equalized PDSCH symbols
hcd = comm.ConstellationDiagram('Title','Equalized PDSCH Symbols',...
                                'ShowReferenceConstellation',false);

%% Transmitter Design: System Architecture
%DL-SCH 下行链路仿真
txsim.RC = 'R.7';       % Base RMC configuration, 10 MHz bandwidth
txsim.NCellID = 88;     % Cell identity
txsim.NFrame = 700;     % Initial frame number
txsim.TotFrames = 1;    % Number of frames to generate
txsim.DesiredCenterFrequency = 2.45e9; % Center frequency in Hz
txsim.NTxAnts = 1;      % Number of transmit antennas

%%
% *Prepare Image File 图像->比特流

fileTx = 'peppers.png';            % Image file name
fData = imread(fileTx);            % Read image data from file
scale = 0.3;                       % Image scaling factor
origSize = size(fData);            % Original input image size
scaledSize = max(floor(scale.*origSize(1:2)),1); % Calculate new image size
heightIx = min(round(((1:scaledSize(1))-0.5)./scale+0.5),origSize(1));
widthIx = min(round(((1:scaledSize(2))-0.5)./scale+0.5),origSize(2));
fData = fData(heightIx,widthIx,:); % Resize image
imsize = size(fData);              % Store new image size
binData = dec2bin(fData(:),8);     % Convert to 8 bit unsigned binary
trData = reshape((binData-'0').',1,[]).'; % Create binary stream

%% 传输图片
% Plot transmit image
figure(1);
 imFig.Visible = 'on';
subplot(211); 
    imshow(fData);
    title('Transmitted Image');
subplot(212);
    title('Received image will appear here...');
    set(gca,'Visible','off'); % Hide axes
    set(findall(gca, 'type', 'text'), 'visible', 'on'); % Unhide title

pause(1); % Pause to plot Tx image
    
%%
% *Generate Baseband LTE Signal*

rmc = lteRMCDL(txsim.RC);

% Calculate the required number of LTE frames based on the size of the
% image data，计算需要多少个LTE帧。
trBlkSize = rmc.PDSCH.TrBlkSizes;
txsim.TotFrames = ceil(numel(trData)/sum(trBlkSize(:)));

% Customize RMC parameters
rmc.NCellID = txsim.NCellID;
rmc.NFrame = txsim.NFrame;
rmc.TotSubframes = txsim.TotFrames*10; % 10 subframes per frame
rmc.CellRefP = txsim.NTxAnts; % Configure number of cell reference ports
rmc.PDSCH.RVSeq = 0;

% Fill subframe 5 with dummy data
rmc.OCNGPDSCHEnable = 'On';
rmc.OCNGPDCCHEnable = 'On';

% If transmitting over two channels enable transmit diversity
% if rmc.CellRefP == 2
%     rmc.PDSCH.TxScheme = 'TxDiversity';
%     rmc.PDSCH.NLayers = 2;
%     rmc.OCNGPDSCH.TxScheme = 'TxDiversity';
% end

fprintf('\nGenerating LTE transmit waveform:\n')
fprintf('  Packing image data into %d frame(s).\n\n', txsim.TotFrames);

% Pack the image data into a single LTE frame
[eNodeBOutput,txGrid,rmc] = lteRMCDLTool(rmc,trData);

% ------------------ 发射机部分结束--------%

% ------------------接收机部分-------------% 
%% Receiver Design: System Architecture
% The general structure of the LTE receiver can be described as follows:
%
% # Capture a suitable number of frames of the transmitted LTE signal using
% SDR hardware.
% # Determine and correct the frequency offset of the received signal.
% # Synchronize the captured signal to the start of an LTE frame.
% # OFDM demodulate the received signal to get an LTE resource grid.
% # Perform a channel estimation for the received signal.
% # Decode the PDSCH and DL-SCH to obtain the transmitted data from the
% transport blocks of each radio frame.
% # Recombine received transport block data to form the received image.
%
% This example plots the power spectral density of the captured waveform,
% and shows visualizations of the estimated channel, equalized PDSCH
% symbols, and received image.
%%
% *Receiver Setup*

% User defined parameters --- configure the same as transmitter
rxsim = struct;
rxsim.RadioFrontEndSampleRate = rmc.SamplingRate; % Configure for same sample rate
                                                       % as transmitter
rxsim.RadioCenterFrequency = txsim.DesiredCenterFrequency;
rxsim.NRxAnts = txsim.NTxAnts;
rxsim.FramesPerBurst = txsim.TotFrames+1; % Number of LTE frames to capture in each burst.
                                          % Capture 1 more LTE frame than transmitted to  
                                          % allow for timing offset wraparound...
rxsim.numBurstCaptures = 1; % Number of bursts to capture

% Derived parameters
samplesPerFrame = 10e-3*rxsim.RadioFrontEndSampleRate; % LTE frames period is 10 ms

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

%%
% *Signal Capture and Processing*
enbDefault = enb;

while rxsim.numBurstCaptures
    % Set default LTE parameters
    enb = enbDefault;
    
%     % SDR Capture
%     fprintf('\nStarting a new RF capture.\n\n')
%     len = 0;
%     for frame = 1:rxsim.FramesPerBurst
%         while len == 0
%             % Store one LTE frame worth of samples
%             [data,len,lostSamples] = step(rx);
%             burstCaptures(:,:,frame) = data;
%         end
%         if lostSamples
%             warning('Dropped samples');
%         end
%         len = 0;
%     end    
%     if rxsim.NRxAnts == 2
%         rxWaveform = reshape(permute(burstCaptures,[1 3 2]), ...
%                         rxsim.FramesPerBurst*samplesPerFrame,rxsim.NRxAnts);
%         hsa.ShowLegend = true; % Turn on legend for spectrum analyzer
%         hsa.ChannelNames = {'SDR Channel 1','SDR Channel 2'};
%     else
%     rxWaveform = burstCaptures(:);
      rxWaveform  = eNodeBOutput;
%     end
    
    % Show power spectral density of captured burst
    hsa.SampleRate = rxsim.RadioFrontEndSampleRate;
    step(hsa,rxWaveform);
    
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

               release(hcd); % Release previous constellation plot
               step(hcd,rxEncodedSymb{:}); % Plot current constellation
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
%        figure(hhest);
%        hhest.Visible = 'On';
%        surf(abs(hest(:,focalFrameIdx,1,1)));
%         shading flat;
%         xlabel('OFDM symbol index'); 
%         ylabel('Subcarrier index');
%         zlabel('Magnitude');   
%         title('Estimate of Channel Magnitude Frequency Repsonse');                
    end
    rxsim.numBurstCaptures = rxsim.numBurstCaptures-1;
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
hBER = comm.ErrorRate;
err = step(hBER, decodedRxDataStream(1:length(trData)), trData);
fprintf('  Bit Error Rate (BER) = %0.5f.\n', err(1));
fprintf('  Number of bit errors = %d.\n', err(2));
fprintf('  Number of transmitted bits = %d.\n',length(trData));

% Recreate image from received data
fprintf('\nConstructing image from received data.\n');
str = reshape(sprintf('%d',decodedRxDataStream(1:length(trData))), 8, []).';
decdata = uint8(bin2dec(str));
receivedImage = reshape(decdata,imsize);

% Plot receive image

figure(1);
subplot(212); 
imshow(receivedImage);
title(sprintf('Received Image: %dx%d Antenna Configuration',txsim.NTxAnts, rxsim.NRxAnts));

%% Things to Try
% By default, the example will use multiple antennas for transmission and
% reception of the LTE waveform. You can modify the transmitter and receiver
% to use a single antenna and decrease the transmitter gain, to 
% observe the difference in the EVM and BER after signal reception and 
% processing. You should also be able to see any errors in the displayed,
% received image.