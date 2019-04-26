%% 802.11ac(TM) Transmit Beamforming
% This example shows how the performance of an IEEE(R) 802.11ac(TM) link
% can be improved by beamforming the transmission when channel state
% information is available at the transmitter.

% Copyright 2015 The MathWorks, Inc.

%% Introduction
% Transmit beamforming focuses energy towards a receiver to improve the SNR
% of a link. In this scheme the transmitter is called a beamformer and the
% receiver is called a beamformee. A steering matrix is used by the
% beamformer to direct the energy to the beamformee. The steering matrix is
% calculated using channel state information obtained through channel
% measurements. In IEEE(R) 802.11ac [ <#23 1> ] these measurements are
% obtained by sounding the channel between beamformer and beamformee. To
% sound the channel the beamformer sends an NDP (Null Data Packet) to the
% beamformee. The beamformee uses the channel information provided by
% sounding to calculate a feedback matrix. This matrix is fedback to the
% beamformer in a compressed format. The beamformer can then use the
% feedback matrix to create a steering matrix and beamform transmissions to
% the beamformee. The process of forming the steering matrix is shown in
% the diagram below.
%
% <<vhtBeamformingFeedback.png>>
%
% In 802.11ac the single user beamformee capability is not mandatory.
% Therefore a multi-antenna transmitter may have to use a different scheme
% to transmit packets to a receiver which cannot act as a beamformee. One
% such scheme is spatial expansion. Spatial expansion allows a number of
% space-time streams to be transmitted on a greater number of transmit
% antennas. Using spatial expansion can provide a small transmit diversity
% gain in channels with flat fading when compared to directly mapping
% space-time streams to transmit antennas [ <#23 2> ].
%
% In this example a 4x2 MIMO configuration is considered between a
% transmitter and receiver, with two space-time streams used for a data
% packet transmission. First the scenario of a receiver which is not
% capable of being a beamformee is considered. A transmission is made using
% spatial expansion and the data symbols are recovered and the signal
% quality measured. To show the benefits of transmit beamforming the data
% packet is then transmitted over the same channel realization, but this
% time using transmit beamforming. The performance of the two schemes are
% then compared. These stages are shown in the diagram below.
%
% <<vhtBeamformingDiagram.png>>

%% Waveform Configuration
% A 4x2 MIMO configuration is used in this example with 2 space-time
% streams.

NumTxAnts = 4;  % Number of transmit antennas
NumSTS = 2;     % Number of space-time streams
NumRxAnts = 2;  % Number of receive antennas

%%
% The format specific configuration of a VHT waveform is described using a
% VHT format configuration object. In this example the waveform is
% configured with a 20 MHz bandwidth and the MIMO configuration specified
% above.

cfgVHT = wlanVHTConfig;
cfgVHT.ChannelBandwidth = 'CBW20';
cfgVHT.APEPLength = 4000;
cfgVHT.NumTransmitAntennas = NumTxAnts;
cfgVHT.NumSpaceTimeStreams = NumSTS;
cfgVHT.MCS = 4; % 16-QAM, rate 3/4

%% Channel Configuration
% In this example a TGac channel model is used with delay profile Model-B.
% The channel realization is controlled with a seed to allow repeatability.

tgac = wlanTGacChannel;
tgac.DelayProfile = 'Model-B';
tgac.ChannelBandwidth = cfgVHT.ChannelBandwidth;
tgac.SampleRate = helperSampleRate(cfgVHT);
tgac.NumReceiveAntennas = NumRxAnts;
tgac.NumTransmitAntennas = NumTxAnts;
tgac.TransmitReceiveDistance = 100; % Meters
tgac.RandomStream = 'mt19937ar with seed';
tgac.Seed = 70; % Seed to allow repeatability

%%
% Noise is added to the time domain waveform at the output of the channel
% with a power, |noisePower|.

noisePower = -37; % dBW

%%
% Setup other objects and variables for simulation.

% Recovery configuration
cfgRec = wlanRecoveryConfig('PilotPhaseTracking','None');

% Indices for extracting fields
ind = wlanFieldIndices(cfgVHT);

% AWGN channel to add noise with a specified noise power. The random
% process controlling noise generation is seeded to allow repeatability.
AWGN = comm.AWGNChannel;
AWGN.RandomStream = 'mt19937ar with seed';
AWGN.Seed = 5;
AWGN.NoiseMethod = 'Variance';
AWGN.Variance = 10^(noisePower/10);

% Calculate the expected noise variance after OFDM demodulation
noiseVar = vhtBeamformingNoiseVariance(noisePower,cfgVHT);

% Number of spatial streams
Nss = NumSTS/(cfgVHT.STBC+1);

% Get the number of occupied subcarriers in VHT fields
[data,pilots] = helperSubcarrierIndices(cfgVHT,'VHT');
Nsd = numel(data);   % Number of data subcarriers
Nsp = numel(pilots); % Number of pilot subcarriers
Nst = Nsd+Nsp; % Total number of occupied subcarriers

% Generate a random PSDU which will be transmitted
rng(0); % Set random state for repeatability
psdu = randi([0 1],cfgVHT.PSDULength*8,1);

%% Transmission with Spatial Expansion
% First a transmission is made using spatial expansion. This type of
% transmission may be made by a multi-antenna transmitter to a receiver
% which is not capable of being a beamformee. The |SpatialMapping| property
% of the format configuration object allows different spatial mapping
% schemes to be selected. In this example the example spatial expansion
% matrix provided in IEEE Std 802.11-2012 Section 2.3.11.1.1.2 [ <#23 3> ]
% is used. Therefore a |'Custom'| spatial mapping is configured. The custom
% spatial mapping matrix is used by assigning the |SpatialMappingMatrix| of
% the format configuration object. This matrix describes the mapping of
% each subcarrier for each space-time stream to all transmit antennas.
% Therefore the size of the spatial mapping matrix used is
% |Nst-by-Nsts-by-Nt|. |Nst| is the number of occupied subcarriers, |Nsts|
% is the number of space-time streams, and |Nt| is the number of transmit
% antennas. The spatial mapping matrix duplicates some of the space-time
% streams to form the desired number of transmit streams.

% Configure a spatial expansion transmission
vhtSE = cfgVHT;
vhtSE.SpatialMapping = 'Custom'; % Use custom spatial expansion matrix
vhtSE.SpatialMappingMatrix = helperSpatialExpansionMatrix(vhtSE);

% Generate waveform
tx = wlanWaveformGenerator(psdu,vhtSE);

% Pass waveform through a fading channel and add noise. Trailing zeros
% are added to allow for channel filter delay.
rx = step(tgac,[tx; zeros(15,NumTxAnts)]);
reset(tgac); % Allow same channel realization to be used subsequently
rx = step(AWGN,rx);
reset(AWGN); % Allow same noise realization to be used subsequently

% Synchronize
tOff = helperSymbolTiming(rx,vhtSE.ChannelBandwidth)-double(ind.LLTF(1))-1;
rx = rx(tOff+1:end,:);

% Channel estimation
vhtltf = rx(ind.VHTLTF(1):ind.VHTLTF(2),:);
vhtltfDemod = wlanVHTLTFDemodulate(vhtltf,vhtSE);
chanEstSE = wlanVHTLTFChannelEstimate(vhtltfDemod,vhtSE);

%%
% The received data field is demodulated and equalized to recover OFDM
% symbols for each spatial stream.

% Demodulate and equalize the data
vhtdata = rx(ind.VHTData(1):ind.VHTData(2),:);
[~,~,symSE] = wlanVHTDataRecover(vhtdata,chanEstSE,noiseVar,vhtSE,cfgRec);

%%
% The constellation of each spatial stream is plotted below.

vhtBeamformingPlotConstellation(symSE, ...
    'Spatial Expansion Transmission Equalized Symbols');

%%
% The variance in the constellation is approximately the same for each
% spatial stream as the SNRs are approximately the same. This is because
% the average power in the channel is on average approximately the same per
% space-time stream:

disp(['Mean received channel power per space-time stream ' ...
    'with spatial expansion: '])
for i = 1:NumSTS
    fprintf('  Space-time stream %d: %2.2f W\n',i, ...
        sum(mean(chanEstSE(:,i,:).*conj(chanEstSE(:,i,:)),1),3))
end

%% Transmission with Beamforming
% When the receiver is capable of being a beamformee, a beamformed
% transmission can create a higher SNR compared to spatial expansion. We
% will now show the advantage of having channel state information available
% to create and use a steering matrix. To calculate a beamforming steering
% matrix, an NDP is passed through the channel. |'Direct'| spatial mapping
% is used for the NDP transmission and the number of space-time streams is
% configured to match the number of transmit antennas. This allows the
% VHT-LTF to be used to sound channels between each of the transmit
% antennas and receive antennas. The calculated beamforming matrix is then
% used to beamform a transmission through the channel. The same channel
% realization is used for sounding and data transmission and there is no
% feedback compression between beamformee and beamformer, therefore the
% beamforming can be regarded as perfect in this example.

% Configure a sounding packet
vhtSound = cfgVHT;
vhtSound.APEPLength = 0; % NDP so no data
vhtSound.NumSpaceTimeStreams = NumTxAnts;
vhtSound.SpatialMapping = 'Direct'; % Each TxAnt carries a STS

% Generate sounding waveform
soundingPSDU = [];
tx = wlanWaveformGenerator(soundingPSDU,vhtSound);

% Pass sounding waveform through the channel and add noise. Trailing zeros
% are added to allow for channel filter delay.
rx = step(tgac,[tx; zeros(15,NumTxAnts)]);
reset(tgac); % Allow same channel realization to be used subsequently
rx = step(AWGN,rx);
reset(AWGN); % Allow same noise realization to be used subsequently

% Synchronize
tOff = helperSymbolTiming(rx,vhtSound.ChannelBandwidth) ...
    -double(ind.LLTF(1))-1;
rx = rx(tOff+1:end,:);

%%
% Channel estimation is performed using the sounding packet to estimate the
% actual channel response between each transmit and receive antenna.

% Channel estimation
vhtLLTFInd = wlanFieldIndices(vhtSound,'VHT-LTF');
vhtltf = rx(vhtLLTFInd(1):vhtLLTFInd(2),:);
vhtltfDemod = wlanVHTLTFDemodulate(vhtltf,vhtSound);
chanEstSound = wlanVHTLTFChannelEstimate(vhtltfDemod,vhtSound);

%%
% The channel estimated using <matlab:doc('wlanVHTLTFChannelEstimate')
% wlanVHTLTFChannelEstimate> includes cyclic shifts applied at the
% transmitter to each space-time stream. To calculate a beamforming
% steering matrix the cyclic shifts applied at the transmitter are removed
% from the channel estimate.

% Remove impact of cyclic shift from channel estimate
chanEstSound = vhtBeamformingRemoveCSD(chanEstSound,vhtSound);

%%
% In this example the beamforming steering matrix is calculated using
% singular value decomposition (SVD) [ <#23 3> ]. The SVD of the channel
% matrix results in two unitary matrices, |U| and |V|, and a diagonal
% matrix of singular values |S|. The first |NumSTS| columns of |V| per
% subcarrier are used as the beamforming steering matrix. The SVD is
% computed using the function <matlab:doc('svd') svd>.

chanEstPerm = permute(chanEstSound,[3 2 1]); % permute to Nr-by-Nt-by-Nst
V = zeros(Nst,NumTxAnts,NumRxAnts);
for i = 1:Nst
    [U,S,V(i,:,:)] = svd(chanEstPerm(:,:,i),'econ');
end
steeringMatrix = V(:,:,1:NumSTS); % Nst-by-Nt-by-Nsts

%%
% The beamforming steering matrix calculated above is applied as a custom
% spatial mapping matrix and is used to send data through the same channel.

% Configure a transmission with beamforming
vhtBF = cfgVHT;
vhtBF.SpatialMapping = 'Custom';
% Permute steering matrix to Nst-by-Nsts-by-Nt
vhtBF.SpatialMappingMatrix = permute(steeringMatrix,[1 3 2]); 

% Generate beamformed data transmission
tx = wlanWaveformGenerator(psdu,vhtBF);

% Pass through the channel and add noise. Trailing zeros
% are added to allow for channel filter delay.
rx = step(tgac,[tx; zeros(15,NumTxAnts)]);
rx = step(AWGN,rx);

% Synchronize
tOff = helperSymbolTiming(rx,vhtBF.ChannelBandwidth)-double(ind.LLTF(1))-1;
rx = rx(tOff+1:end,:);

% Channel estimation
vhtltf = rx(ind.VHTLTF(1):ind.VHTLTF(2),:);
vhtltfDemod = wlanVHTLTFDemodulate(vhtltf,vhtBF);
chanEstBF = wlanVHTLTFChannelEstimate(vhtltfDemod,vhtBF);

%%
% The received data field is demodulated and equalized to recover OFDM
% symbols for each spatial stream.

% Demodulate and equalize the data
vhtdata = rx(ind.VHTData(1):ind.VHTData(2),:);
[~,~,symBF] = wlanVHTDataRecover(vhtdata,chanEstBF,noiseVar,vhtBF,cfgRec);

%%
% The equalized constellation for each spatial stream is plotted below.
% Note that the higher order spatial stream has a larger variance. This is
% due to the ordered singular values of the channels used in SVD
% beamforming.

vhtBeamformingPlotConstellation(symBF, ...
    'Beamformed Transmission Equalized Symbols');

%%
% This ordering is also visible in the average power of the received
% space-time streams. The power of the received first space-time stream is
% larger than the second space-time stream. This is because the received
% signal strength is a function of the singular values of the channel which
% SVD orders in a decreasing fashion.

disp(['Mean received channel power per space-time stream ' ...
    'with SVD transmit beamforming: '])
for i = 1:NumSTS
    fprintf('  Space-time stream %d: %2.2f W\n',i, ...
        sum(mean(chanEstBF(:,i,:).*conj(chanEstBF(:,i,:)),1),3))
end

%% Comparison and Conclusion
% The figure below plots the equalized constellation from the spatial
% expansion and beamformed transmissions for all spatial streams. Note the
% improved constellation using SVD-based transmit beamforming.

figure
plot(symSE(:),'Marker','.','LineStyle','none','Color',[0.929 0.694 0.125])
hold on
plot(symBF(:),'Marker','.','LineStyle','none','Color',[0.494 0.1840 0.556])
title('Equalized Symbol Comparison')
str = sprintf('%dx%d',NumTxAnts,NumRxAnts);
legend([str ' Spatial Expansion'],[str ' Transmit Beamforming'])
xlabel('Real')
ylabel('Imag')
xlim([-2 2])
ylim([-2 2])

%% 
% The improvement can also be measured through the RMS and maximum error
% vector magnitude (EVM). EVM is a measure of demodulated signal quality.

EVM = comm.EVM;
EVM.AveragingDimensions = [1 2]; % Average over all subcarriers and symbols
EVM.MaximumEVMOutputPort = true;

refSE = helperClosestConstellationPoint(symSE,vhtSE);
[rmsEVMSE,maxEVMSE] = step(EVM,refSE,symSE); % EVM using spatial expansion
refBF = helperClosestConstellationPoint(symBF,vhtBF);
[rmsEVMBF,maxEVMBF] = step(EVM,refBF,symBF); % EVM using beamforming

for i = 1:Nss
    fprintf(['Spatial stream %d EVM:\n' ...
        '  Spatial expansion:    %2.1f%% RMS, %2.1f%% max\n' ...
        '  Transmit beamforming: %2.1f%% RMS, %2.1f%% max\n'], ...
        i,rmsEVMSE(i),maxEVMSE(i),rmsEVMBF(i),maxEVMBF(i));
end

%%
% This example demonstrates that if a receiver is capable of being a
% beamformee, the SNR can potentially be improved when a transmission is
% beamformed compared to a spatial expansion transmission. The increase in
% received power when using beamforming can lead to more reliable
% demodulation or potentially even a higher order modulation and coding
% scheme to be used for the transmission.
%
% In a realistic operational simulation the performance of beamforming
% would be degraded due to the delay between channel state information
% calculation and feedback by the beamformee and feedback quantization. For
% more information see  [ <#23 2> ].

%% Appendix
% This example uses the following helper functions:
%
% * <matlab:edit('helperSampleRate.m') helperSampleRate.m>
% * <matlab:edit('helperSubcarrierIndices.m') helperSubcarrierIndices.m>
% * <matlab:edit('vhtBeamformingNoiseVariance.m') vhtBeamformingNoiseVariance.m>
% * <matlab:edit('helperSpatialExpansionMatrix.m') helperSpatialExpansionMatrix.m>
% * <matlab:edit('vhtBeamformingPlotConstellation.m') vhtBeamformingPlotConstellation.m>
% * <matlab:edit('vhtBeamformingRemoveCSD.m') vhtBeamformingRemoveCSD.m>
% * <matlab:edit('helperClosestConstellationPoint.m') helperClosestConstellationPoint.m>

%% Selected Bibliography
% # IEEE Std 802.11ac(TM)-2013 IEEE Standard for Information technology -
% Telecommunications and information exchange between systems - Local and
% metropolitan area networks - Specific requirements - Part 11: Wireless
% LAN Medium Access Control (MAC) and Physical Layer (PHY) Specifications -
% Amendment 4: Enhancements for Very High Throughput for Operation in Bands
% below 6 GHz.
% # Perahia, Eldad, and Robert Stacey. Next Generation Wireless LANS:
% 802.11n and 802.11ac. Cambridge University Press, 2013.
% # IEEE Std 802.11(TM)-2012 IEEE Standard for Information technology -
% Telecommunications and information exchange between systems - Local and
% metropolitan area networks - Specific requirements - Part 11: Wireless
% LAN Medium Access Control (MAC) and Physical Layer (PHY) Specifications.


displayEndOfDemoMessage(mfilename)