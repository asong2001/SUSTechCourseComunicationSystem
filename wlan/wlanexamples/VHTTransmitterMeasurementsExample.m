%% 802.11ac(TM) Transmitter Modulation Accuracy and Spectral Emission Testing
%
% This example shows how to perform transmitter modulation accuracy and
% spectrum emission mask and flatness measurements on an IEEE(R)
% 802.11ac(TM) waveform.

% Copyright 2015 The MathWorks, Inc.

%% Introduction
% The transmitter modulation accuracy, required spectrum mask and required
% spectral flatness for given configurations are specified in Section
% 22.3.18 of the 802.11ac standard [ <#18 1> ]. This example shows how
% these measurements can be performed on a waveform. The waveform is
% generated with WLAN System Toolbox(TM), but a waveform captured with a
% spectrum analyzer could be used.
%
% A waveform consisting of five 80 MHz VHT packets separated by 10
% microsecond gaps is generated. Random data is used in each packet and
% 256QAM modulation is used. The baseband waveform is upsampled and
% filtered to reduce the out of band emissions to meet the spectral mask
% requirement. A high power amplifier (HPA) model is used, which introduces
% inband distortion and spectral regrowth. The spectral emission mask
% measurement is performed on the upsampled waveform after the high power
% amplifier modeling. The waveform is downsampled and the error vector
% magnitude (EVM) of the VHT Data field is measured to determine the
% modulation accuracy. The spectral flatness is additionally measured. The
% example is illustrated in the following diagram:
%
% <<VHTTransmitterMeasurementsDiagram.png>>

%% IEEE 802.11ac VHT Packet Configuration
% In this example an IEEE 802.11ac waveform consisting of multiple VHT
% format packets is generated. The format specific configuration of a VHT
% waveform is described using a VHT format configuration object. The object
% is created using the <matlab:doc('wlanVHTConfig') wlanVHTConfig>
% function. The properties of the object contain the configuration. In this
% example the object is configured for a 80 MHz bandwidth. One spatial
% stream is transmitted per antenna to allow the modulation accuracy to be
% measured per spatial stream, therefore no space time block coding is
% used.

cfgVHT = wlanVHTConfig;            % Create packet configuration
cfgVHT.ChannelBandwidth = 'CBW80'; % 80 MHz
cfgVHT.NumTransmitAntennas = 1;    % One transmit antenna
cfgVHT.NumSpaceTimeStreams = 1;    % One space-time stream
cfgVHT.STBC = false;               % No STBC so one spatial stream
cfgVHT.MCS = 8;                    % Modulation: 256 QAM
cfgVHT.APEPLength = 3000;          % A-MPDU length pre-EOF padding in bytes

%% Baseband Waveform Generation
% The waveform generator <matlab:doc('wlanWaveformGenerator')
% wlanWaveformGenerator> can be configured to generate one or more packets
% and add an idle time between each packet. In this example five packets
% with a 10 microsecond idle period will be created.

numPackets = 5;   % Generate 5 packets
idleTime = 10e-6; % 10 microsecond idle time between packets

%% 
% Random bits for all packets |data| are created and passed as an argument
% to <matlab:doc('wlanWaveformGenerator') wlanWaveformGenerator> along with
% the VHT packet configuration object |cfgVHT|. This configures the
% waveform generator to synthesize an 802.11ac VHT waveform. The waveform
% generator is additionally configured using name-value pairs to generate
% multiple packets with a specified idle time between each packet.

% Create random data; PSDULength is in bytes
savedState = rng(0); % Set random state
data = randi([0 1],cfgVHT.PSDULength*8*numPackets,1);

% Generate a multi-packet waveform
txWaveform = wlanWaveformGenerator(data,cfgVHT, ...
    'NumPackets',numPackets,'IdleTime',idleTime);

% Get the sampling rate of the waveform
fs = helperSampleRate(cfgVHT);
disp(['Baseband sampling rate: ' num2str(fs/1e6) ' Msps']);

%% Oversampling and Filtering
% Spectral filtering is used to reduce the out of band spectral emissions
% owing to the implicit rectangular pulse shaping in the OFDM modulation,
% and spectral regrowth caused by the high power amplifier model. To model
% the effect of a high power amplifier on the waveform and view the out of
% band spectral emissions the waveform must be oversampled. Oversampling
% requires an interpolation filter to remove spectral images caused by
% upsampling. In this example the waveform is oversampled with an
% interpolation filter which also acts as a spectral filter. This allows
% the waveform to meet spectral mask requirements. The waveform is
% oversampled and filtered using <matlab:doc('dsp.FIRInterpolator')
% dsp.FIRInterpolator>.

% Oversample the waveform
osf = 3;         % Oversampling factor
filterLen = 120; % Filter length
beta = 0.5;      % Design parameter for Kaiser window

% Generate filter coefficients
coeffs = osf.*firnyquist(filterLen,osf,kaiser(filterLen+1,beta)); 
coeffs = coeffs(1:end-1); % Remove trailing zero
FIRINTERP = dsp.FIRInterpolator(osf,'Numerator',coeffs);
txWaveform = step(FIRINTERP,txWaveform);

% Plot the magnitude and phase response of the filter applied after
% oversampling
h = fvtool(FIRINTERP);
h.Analysis = 'freq';           % Plot magnitude and phase responses
h.FS = osf*fs;                 % Set sampling rate
h.NormalizedFrequency = 'off'; % Plot responses against frequency

%% High Power Amplifier Modeling
% The high power amplifier introduces nonlinear behavior in the form of
% inband distortion and spectral regrowth. The Rapp model is used to
% simulate power amplifiers in 802.11ac [ <#18 2> ]. The Rapp model causes
% AM/AM distortion and is modeled with
% <matlab:doc('comm.MemorylessNonlinearity') comm.MemorylessNonlinearity>.
% The high power amplifier is backed-off to operate below the saturation
% point to reduce distortion. The backoff is controlled by the variable
% |hpaBackoff|.

hpaBackoff = 8; % dB

% Create and configure a memoryless nonlinearity to model the amplifier
HPA = comm.MemorylessNonlinearity;
HPA.Method = 'Rapp model';
HPA.Smoothness = 3; % p parameter
HPA.LinearGain = -hpaBackoff;

% Apply the model to each transmit antenna
for i=1:cfgVHT.NumTransmitAntennas
    txWaveform(:,i) = step(HPA,txWaveform(:,i));
end

%% 
% Thermal noise is added to the waveform with a 6 dB noise figure [ <#18 3> ].

NF = 6;         % Noise figure (dB)
BW = fs*osf;    % Bandwidth (Hz)
k = 1.3806e-23; % Boltzman constant (J/K)
T = 290;        % Ambient temperature (K)
noisePower = 10*log10(k*T*BW)+NF;

AWGN = comm.AWGNChannel('NoiseMethod','Variance', ...
    'Variance',10^(noisePower/10));
txWaveform = step(AWGN,txWaveform);

%% Modulation Accuracy (EVM) and Spectral Flatness Measurements
% The oversampled waveform is resampled down to baseband for physical layer
% processing and EVM and spectral flatness measurements. As part of the
% resampling a low-pass anti-aliasing filter is applied before
% downsampling. The impact of the low-pass filter will be visible in the
% spectral flatness measurement. The waveform is resampled to baseband
% using <matlab:doc('dsp.FIRDecimator') dsp.FIRDecimator> with the same
% coefficients used for oversampling earlier in the example.

% Resample to baseband
FIRDEC = dsp.FIRDecimator(osf,'Numerator',coeffs);
rxWaveform = step(FIRDEC,txWaveform);

%%
% Each packet within |rxWaveform| is detected, synchronized and extracted.
% The EVM and spectral flatness measurements are made for each packet. The
% following steps are performed for each packet:
%
% * The start of the packet is detected
% * The non-HT fields are extracted and coarse carrier frequency offset
% (CFO) estimation and correction are performed
% * The frequency corrected non-HT fields are used to estimate fine symbol
% timing
% * The packet is extracted from the waveform using the fine symbol timing
% offset
% * The extracted packet is corrected with the coarse CFO estimate
% * The L-LTF is extracted and used to estimate the fine CFO. The offset is
% corrected for the whole packet
% * The L-LTF is used to estimate the noise power
% * The VHT-LTF is extracted and channel estimation is performed for each
% of the transmit streams
% * The channel estimate is used to measure the spectral flatness
% * The VHT Data field is extracted, OFDM demodulated, phase corrected and
% equalized using the channel estimate
% * For each data-carrying subcarrier in each spatial stream, the closest
% constellation point is found and the EVM computed
%
% The processing chain is shown in the diagram below:
%
% <<VHTTransmitterMeasurementsRxChain.png>>
%
% Note the VHT-LTF symbols include pilot symbols to allow for phase
% tracking, but this is not done in this example.
%
% The spectral flatness is tested for each packet by measuring the
% deviation in the magnitude of individual subcarriers in the channel
% estimate against the average [ <#18 1> ]. These deviations are plotted
% for each packet using the helper function
% |vhtTxSpectralFlatnessMeasurement|. The average EVM per data-carrying
% subcarrier, and the equalized symbols are plotted for each packet.
%
% The function <matlab:doc('wlanVHTDataRecover') wlanVHTDataRecover> is
% used to demodulate, equalize and decode the VHT Data symbols. The
% equalized symbols are used in this example to measure the modulation
% accuracy. This function is parameterized using a
% <matlab:doc('wlanRecoveryConfig') wlanRecoveryConfig> object. The object
% is parameterized to perform pilot phase tracking and zero forcing
% equalization as required by the standard.

% Configure VHT Data symbol recovery
cfgRec = wlanRecoveryConfig;
cfgRec.EqualizationMethod = 'ZF';    % Use zero forcing algorithm
cfgRec.PilotPhaseTracking = 'PreEQ'; % Use pilot phase tracking

%%
% Two different EVM measurements are made in this example using two
% instances of <matlab:doc('comm.EVM') comm.EVM>. The first measurement is
% the RMS EVM per packet. For this measurement the EVM is averaged over
% subcarriers, OFDM symbols and spatial streams.

EVMPerPkt = comm.EVM;
EVMPerPkt.AveragingDimensions = [1 2 3]; % Nst-by-Nsym-by-Nss
EVMPerPkt.Normalization = 'Average constellation power';

%%
% The second measurement is the RMS EVM per subcarrier per spatial stream
% for a packet. As spatial streams are mapped directly to antennas in this
% setup, this measurement can help detect frequency dependent impairments
% which may affect individual RF chains differently. For this measurement
% the EVM is only averaged over OFDM symbols.

% Measure average EVM over symbols
EVMPerSC = comm.EVM;
EVMPerSC.AveragingDimensions = 2; % Nst-by-Nsym-by-Nss
EVMPerSC.Normalization = 'Average constellation power';

%%
% The following code configures objects and variables for processing.

% Indices for accessing each field within the time-domain packet
ind = wlanFieldIndices(cfgVHT);

% Create and configure an object to correct the estimated carrier frequency
% offset in the waveform.
PFO = comm.PhaseFrequencyOffset;
PFO.SampleRate = fs;
PFO.PhaseOffset = 0;
PFO.FrequencyOffsetSource = 'Input port';

rxWaveformLength = size(rxWaveform,1);
pktLength = double(ind.VHTData(2));

% Minimum length of data we can detect; length of the L-STF in samples
minPktLen = double(ind.LSTF(2)-ind.LSTF(1))+1;

% Setup the measurement plots
[hSF,hCon,hEVM] = vhtTxSetupPlots(cfgVHT);

rmsEVM = [];         % Empty as unknown number of packets
pktOffsetStore = []; % Empty as unknown number of packets

rng(savedState); % Restore random state

%%
% A <matlab:doc('while') while> loop is used to detect and process packets
% within the received waveform. The sample offset |searchOffset| is used to
% index into |rxWaveform| to detect a packet. The first packet within
% |rxWaveform| is detected and processed. The sample index offset
% |searchOffset| is then incremented to move beyond the processed packet in
% |rxWaveform| and the next packet is detected and processed until no
% further packets are detected.

pktNum = 0;
searchOffset = 0; % Start at first sample (no offset)
while (searchOffset+minPktLen)<=rxWaveformLength
    % Packet detect
    pktStartIdx = helperPacketDetect(rxWaveform(1+searchOffset:end,:), ...
        cfgVHT.ChannelBandwidth);
    % Packet offset from start of waveform
    pktOffset = searchOffset+pktStartIdx-1; 
    % If no packet detected or offset outwith bounds of waveform then stop
    if isempty(pktStartIdx) || (pktOffset<0) || ...
            ((pktOffset+ind.LSIG(2))>rxWaveformLength)
        break;
    end
    
    % Extract non-HT fields and perform coarse frequency offset correction
    % to allow for reliable symbol timing
    nonht = rxWaveform(pktOffset+(ind.LSTF(1):ind.LSIG(2)),:);  
    coarsefreqOff = wlanCoarseCFOEstimate(nonht,cfgVHT.ChannelBandwidth);
    nonht = step(PFO,nonht,-coarsefreqOff);
    release(PFO); % Allow for subsequent processing with a different size
    
    % Determine start of L-LTF
    lltfIdx = helperSymbolTiming(nonht,cfgVHT.ChannelBandwidth);
    % If no L-LTF detected skip samples and continue searching
    if isempty(lltfIdx)
        searchOffset = pktOffset+double(ind.LSTF(2))+1;
        continue; 
    end
    % Determine packet offset given the offset between the expected start
    % of L-LTF and actual start of L-LTF
    pktOffset = pktOffset+lltfIdx-double(ind.LLTF(1));
    % If offset is without bounds of waveform  skip samples and continue
    % searching within remainder of the waveform
    if (pktOffset<0) || ((pktOffset+pktLength)>rxWaveformLength)
        searchOffset = pktOffset+double(ind.LSTF(2))+1;
        continue;
    end  
    
    % Timing synchronization complete; extract the detected packet
    rxPacket = rxWaveform(pktOffset+(1:pktLength),:);
    pktNum = pktNum+1;
    disp(['  Packet ' num2str(pktNum) ' at index: ' num2str(pktOffset+1)]);
    
    % Apply coarse frequency correction to the extracted packet
    rxPacket = step(PFO,rxPacket,-coarsefreqOff);
    release(PFO); % Allow for subsequent processing with a different size
    
    % Perform fine frequency offset correction on the extracted packet
    lltf = rxPacket(ind.LLTF(1):ind.LLTF(2),:); % Extract L-LTF
    fineFreqOff = wlanFineCFOEstimate(lltf,cfgVHT.ChannelBandwidth);
    rxPacket = step(PFO,rxPacket,-fineFreqOff);
    release(PFO); % Release object for subsequent processing
    
    % Estimate noise power in VHT fields
    lltf = rxPacket(ind.LLTF(1):ind.LLTF(2),:);
    demodLLTF = wlanLLTFDemodulate(lltf,cfgVHT);
    noiseVarVHT = helperNoiseEstimate(demodLLTF,cfgVHT);
    
    % Extract VHT-LTF samples, demodulate and perform channel estimation
    vhtltf = rxPacket(ind.VHTLTF(1):ind.VHTLTF(2),:);
    vhtltfDemod = wlanVHTLTFDemodulate(vhtltf,cfgVHT);
    chanEst = wlanVHTLTFChannelEstimate(vhtltfDemod,cfgVHT);
    
    % Spectral flatness measurement
    vhtTxSpectralFlatnessMeasurement(chanEst,cfgVHT,pktNum,hSF);
    
    % Extract VHT Data samples and perform OFDM demodulation, equalization
    % and phase tracking
    vhtdata = rxPacket(ind.VHTData(1):ind.VHTData(2),:);
    [~,~,eqSym] = wlanVHTDataRecover(vhtdata,chanEst,noiseVarVHT, ...
                    cfgVHT,cfgRec);

    % Find the closest constellation point to each equalized symbol
    refSym = helperClosestConstellationPoint(eqSym,cfgVHT);

    % Compute RMS EVM over all spatial streams for packet
    rmsEVM(pktNum) = step(EVMPerPkt,refSym,eqSym); %#ok<SAGROW>
    fprintf('    RMS EVM: %2.2f%%, %2.2fdB\n', ...
        rmsEVM(pktNum),20*log10(rmsEVM(pktNum)/100));

    % Compute RMS EVM per subcarrier and spatial stream for the packet
    evmPerSC = step(EVMPerSC,refSym,eqSym); % Nst-by-1-by-Nss

    % Plot RMS EVM per subcarrier and equalized constellation
    vhtTxEVMConstellationPlots(eqSym,evmPerSC,cfgVHT,pktNum,hCon,hEVM);
    
    % Store the offset of each packet within the waveform
    pktOffsetStore = [pktOffsetStore; pktOffset]; %#ok<AGROW>
    
    % Increment waveform offset and search remaining waveform for a packet
    searchOffset = pktOffset+pktLength+minPktLen;
end

if pktNum>0
fprintf('Average EVM for %d packets: %2.2f%%, %2.2fdB\n', ...
    pktNum,mean(rmsEVM),20*log10(mean(rmsEVM)/100));
else
    disp('No complete packet detected');
end

%% Transmit Spectrum Emission Mask Measurement
% In this example the spectrum emission mask of the filtered and impaired
% waveform after high power amplifier modeling is measured.
%
% A time gated spectral measurement of the VHT Data field is used for the
% transmitter spectrum emission mask test [ <#18 4> ]. As part of the
% baseband processing the start index of each packet within the baseband
% waveform was stored. These indices are used to extract the VHT Data field
% of each packet from the oversampled |txWaveform|. Any delay introduced in
% the baseband processing chain used to determine the packet indices must
% be accounted for when gating the VHT data field within |txWaveform|. The
% extracted VHT Data fields are concatenated in preparation for
% measurement.

startIdx = osf*(ind.VHTData(1)-1)+1; % Upsampled start of VHT Data
endIdx = osf*ind.VHTData(2);         % Upsampled end of VHT Data
delay = grpdelay(FIRDEC,1);          % Group delay of downsampling filter
idx = zeros(endIdx-startIdx+1,pktNum);
for i=1:pktNum
    % Start of packet in txWaveform
    pktOffset = osf*pktOffsetStore(i)-delay;
    % Indices of VHT Data in txWaveform
    idx(:,i) = (pktOffset+(startIdx:endIdx));
end
gatedVHTData = txWaveform(idx(:),:);

%%
% The spectral mask is specified by the standard relative to the peak power
% spectral density. The plot generated by the helper function
% <matlab:edit('helperSpectralMaskTest.m') helperSpectralMaskTest> overlays
% the required mask with the measured PSD.

helperSpectralMaskTest(gatedVHTData,fs,osf);

%% Conclusion and Further Exploration
% Four results are plotted by this example; spectral flatness, RMS EVM per
% subcarrier, equalized constellation, and spectral mask.
%
% The high power amplifier model introduces significant inband distortion
% and spectral regrowth which is visible in the EVM results, noisy
% constellation and out-of-band emissions in the spectral mask plot. Try
% increasing the high power amplifier backoff and note the improved EVM,
% constellation and lower out-of-band emissions.
% 
% The spectral filtering and downsampling (to bring the waveform to
% baseband for processing) stages include filtering. These filter responses
% affect the spectral flatness measurement. The ripple in the spectral
% flatness measurement is mainly due to downsampling to baseband. Try using
% different filters or filter lengths and note the impact on the spectral
% flatness.

%% Appendix
% This example uses the following helper functions:
%
% * <matlab:edit('vhtTxSetupPlots.m') vhtTxSetupPlots.m>
% * <matlab:edit('vhtTxSpectralFlatnessMeasurement.m') vhtTxSpectralFlatnessMeasurement.m>
% * <matlab:edit('helperClosestConstellationPoint.m') helperClosestConstellationPoint.m>
% * <matlab:edit('vhtTxEVMConstellationPlots.m') vhtTxEVMConstellationPlots.m>
% * <matlab:edit('helperSpectralMaskTest.m') helperSpectralMaskTest.m>
% * <matlab:edit('helperPacketDetect.m') helperPacketDetect.m>
% * <matlab:edit('helperSymbolTiming.m') helperSymbolTiming.m>
% * <matlab:edit('helperSampleRate.m') helperSampleRate.m>
% * <matlab:edit('helperNoiseEstimate.m') helperNoiseEstimate.m>

%% Selected Bibliography
% # IEEE Std 802.11ac(TM)-2013 IEEE Standard for Information technology -
% Telecommunications and information exchange between systems - Local and
% metropolitan area networks - Specific requirements - Part 11: Wireless
% LAN Medium Access Control (MAC) and Physical Layer (PHY) Specifications -
% Amendment 4: Enhancements for Very High Throughput for Operation in Bands
% below 6 GHz.
% # Loc and Cheong. IEEE P802.11 Wireless LANs. TGac Functional
% Requirements and Evaluation Methodology Rev. 16. 2011-01-19.
% # Perahia, E., and R. Stacey. Next Generation Wireless LANs: 802.11n and
% 802.11ac. 2nd Edition. United Kingdom: Cambridge University Press, 2013.
% # Archambault, Jerry, and Shravan Surineni. "IEEE 802.11 spectral
% measurements using vector signal analyzers." RF Design 27.6 (2004):
% 38-49.

displayEndOfDemoMessage(mfilename)