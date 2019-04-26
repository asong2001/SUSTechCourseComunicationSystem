%% 802.11ac(TM) Packet Error Rate Simulation for 8x8 TGac Channel
%
% This example shows how to measure the packet error rate of an IEEE(R)
% 802.11ac(TM) VHT link using an end-to-end simulation with a fading TGac
% channel model and additive white Gaussian noise.

% Copyright 2015 The MathWorks, Inc.

%% Introduction
% In this example an end-to-end simulation is used to determine the packet
% error rate for an 802.11ac [ <#12 1> ] VHT link with a fading channel at
% a selection of SNR points. At each SNR point multiple packets are
% transmitted through a channel, demodulated and the PSDUs recovered. The
% PSDUs are compared to those transmitted to determine the number of packet
% errors and hence the packet error rate. Packet detection, timing
% synchronization, carrier frequency offset correction and phase tracking
% are performed by the receiver. The processing for each packet is
% summarized in the following
%
% <<VHTMIMOPERDiagram.png>>
%
% This example also demonstrates how a <matlab:doc('parfor') parfor> loop
% can be used instead of the <matlab:doc('for') for> loop when simulating
% each SNR point to speed up a simulation. <matlab:doc('parfor') parfor>,
% as part of the Parallel Computing Toolbox(TM), executes processing for
% each SNR in parallel to reduce the total simulation time.

%% Waveform Configuration
% An 802.11ac VHT transmission is simulated in this example. A VHT format
% configuration object contains the format specific configuration of the
% transmission. The object is created using the
% <matlab:doc('wlanVHTConfig') wlanVHTConfig> function. The properties of
% the object contain the configuration. In this example the object is
% configured for a 80 MHz channel bandwidth, 8 transmit antennas, 8
% space-time streams, no space time block coding and 256-QAM rate-5/6 (MCS
% 9).

% Create a format configuration object for a 8-by-8 VHT transmission
cfgVHT = wlanVHTConfig;
cfgVHT.ChannelBandwidth = 'CBW80'; % 80 MHz channel bandwidth
cfgVHT.NumTransmitAntennas = 8;    % 8 transmit antennas
cfgVHT.NumSpaceTimeStreams = 8;    % 8 space-time streams
cfgVHT.APEPLength = 3000;          % APEP length in bytes
cfgVHT.MCS = 9;                    % 256-QAM rate-5/6

%% Channel Configuration
% In this example a TGac N-LOS channel model is used with delay profile
% Model-D. For Model-D when the distance between transmitter and receiver
% is greater than or equal to 10 meters, the model is NLOS. This is
% described further in <matlab:doc('wlanTGacChannel') wlanTGacChannel>. An
% 8x8 MIMO channel is simulated in this example therefore 8 receive
% antennas are specified.

% Create and configure the channel
TGac = wlanTGacChannel;
TGac.DelayProfile = 'Model-D';
TGac.NumReceiveAntennas = 8;
TGac.TransmitReceiveDistance = 10; % Distance in meters for NLOS operation
TGac.ChannelBandwidth = cfgVHT.ChannelBandwidth;
TGac.NumTransmitAntennas = cfgVHT.NumTransmitAntennas;
TGac.LargeScaleFadingEffect = 'None';

%% Simulation Parameters
% For each SNR point in the vector |snr| a number of packets are
% generated, passed through a channel and demodulated to determine the
% packet error rate.

snr = 40:5:50;

%%
% The number of packets tested at each SNR point is controlled by two
% parameters:
%
% # |maxNumErrors| is the maximum number of packet errors simulated at each
% SNR point. When the number of packet errors reaches this limit, the
% simulation at this SNR point is complete.
% # |maxNumPackets| is the maximum number of packets simulated at each SNR
% point and limits the length of the simulation if the packet error limit
% is not reached. 
%
% The numbers chosen in this example will lead to a very short simulation.
% For meaningful results we recommend increasing the numbers.

maxNumErrors = 10;   % The maximum number of packet errors at an SNR point
maxNumPackets = 100; % Maximum number of packets at an SNR point

%% Simulation Setup
% In this simulation frequency offset estimation is used.
% <matlab:doc('comm.PhaseFrequencyOffset') comm.PhaseFrequencyOffset>
%  is used to correct the carrier frequency offset estimated
% by <matlab:doc('wlanCoarseCFOEstimate') wlanCoarseCFOEstimate> and
% <matlab:doc('wlanFineCFOEstimate') wlanFineCFOEstimate>.

% Get the baseband sampling rate
fs = helperSampleRate(cfgVHT);

% Phase offset compensation object
PFO = comm.PhaseFrequencyOffset;
PFO.FrequencyOffsetSource = 'Input port';
PFO.SampleRate = fs;

%%
% Set the remaining parameters for the simulation.

% Get the number of occupied subcarriers in VHT fields and FFT length
[vhtData,vhtPilots] = helperSubcarrierIndices(cfgVHT,'VHT');
Nst_vht = numel(vhtData)+numel(vhtPilots);
Nfft = helperFFTLength(cfgVHT);     % FFT length

% Set the sampling rate of the channel
TGac.SampleRate = fs;

% Indices for accessing each field within the time-domain packet
ind = wlanFieldIndices(cfgVHT);

rng(1); % Set random state for repeatability

%% Processing SNR Points
% For each SNR point a number of packets are tested and the packet error
% rate calculated.
%
% For each packet the following processing steps occur:
%
% # A PSDU is created and encoded to create a single packet waveform.
% # The waveform is passed through a different realization of the TGac
% channel model.
% # AWGN is added to the received waveform to create the desired average
% SNR per subcarrier after OFDM demodulation.
% <matlab:doc('comm.AWGNChannel') comm.AWGNChannel> is configured to
% provide the correct SNR. The configuration accounts for normalization
% within the channel by the number of receive antennas, and the noise
% energy in unused subcarriers which are removed during OFDM demodulation.
% # The packet is detected.
% # Coarse carrier frequency offset is estimated and corrected.
% # Fine timing synchronization is established. The L-STF, L-LTF and L-SIG
% samples are provided for fine timing to allow for packet detection at the
% start or end of the L-STF.
% # Fine carrier frequency offset is estimated and corrected.
% # The VHT-LTF is extracted from the synchronized received waveform. The
% VHT-LTF is OFDM demodulated and channel estimation is performed.
% # The VHT Data field is extracted from the synchronized received
% waveform. The PSDU is recovered using the extracted field and the channel
% estimate.
%
% A <matlab:doc('parfor') parfor> loop can be used to parallelize
% processing of the SNR points, therefore for each SNR point an AWGN
% channel is created and configured with <matlab:doc('comm.AWGNChannel')
% comm.AWGNChannel>. To enable the use of parallel computing for increased
% speed comment out the 'for' statement and uncomment the 'parfor'
% statement below.

S = numel(snr);
packetErrorRate = zeros(S,1);
%parfor i = 1:S % Use 'parfor' to speed up the simulation
for i = 1:S     % Use 'for' to debug the simulation
    % Create an instance of the AWGN channel per SNR point simulated
    AWGN = comm.AWGNChannel;
    AWGN.NoiseMethod = 'Signal to noise ratio (SNR)';
    AWGN.SignalPower = 1/TGac.NumReceiveAntennas; % Normalization
    AWGN.SNR = snr(i)-10*log10(Nfft/Nst_vht); % Account for energy in nulls

    % Loop to simulate multiple packets
    numPacketErrors = 0;
    numPkt = 1; % Index of packet transmitted
    while numPacketErrors<=maxNumErrors && numPkt<=maxNumPackets
        % Generate a packet waveform
        txPSDU = randi([0 1],cfgVHT.PSDULength*8,1); % PSDULength in bytes
        tx = wlanWaveformGenerator(txPSDU,cfgVHT);
        
        % Add trailing zeros to allow for channel delay
        tx = [tx; zeros(50,cfgVHT.NumTransmitAntennas)]; %#ok<AGROW>

        % Pass through TGac fading channel model
        rx = step(TGac, tx);
        reset(TGac); % Reset channel to create a different realization

        % Add noise
        rx = step(AWGN,rx);

        % Packet detect
        pktStartIdx = helperPacketDetect(rx,cfgVHT.ChannelBandwidth);
        if isempty(pktStartIdx) % If empty no L-STF detected; packet error
            numPacketErrors = numPacketErrors+1;
            numPkt = numPkt+1;
            continue; % Go to next loop iteration
        end
        pktOffset = pktStartIdx-1; % Packet offset from start of waveform
        
        % Extract L-STF and perform coarse frequency offset correction
        lstf = rx(pktOffset+(ind.LSTF(1):ind.LSTF(2)),:); 
        coarseFreqOff = wlanCoarseCFOEstimate(lstf,cfgVHT.ChannelBandwidth);
        rx = step(PFO,rx,-coarseFreqOff);
        release(PFO); % Release object for subsequent processing
        
        % Extract the Non-HT fields and determine start of L-LTF
        nonhtfields = rx(pktOffset+(ind.LSTF(1):ind.LSIG(2)),:); 
        lltfIdx = helperSymbolTiming(nonhtfields,cfgVHT.ChannelBandwidth);
        
        % Synchronize the received waveform given the offset between the
        % expected start of the L-LTF and actual start of L-LTF
        pktOffset = pktOffset+lltfIdx-double(ind.LLTF(1));
        % If no L-LTF detected or if packet detected outwith the range of
        % expected delays from the channel modeling; packet error
        if isempty(lltfIdx) || pktOffset<0 || pktOffset>50
            numPacketErrors = numPacketErrors+1;
            numPkt = numPkt+1;
            continue; % Go to next loop iteration
        end        
        rx = rx(1+pktOffset:end,:);
             
        % Extract L-LTF and perform fine frequency offset correction
        lltf = rx(ind.LLTF(1):ind.LLTF(2),:); 
        fineFreqOff = wlanFineCFOEstimate(lltf,cfgVHT.ChannelBandwidth);
        rx = step(PFO,rx,-fineFreqOff);
        release(PFO); % Release object for subsequent processing
        
        % Estimate noise power in VHT fields
        lltf = rx(ind.LLTF(1):ind.LLTF(2),:); 
        demodLLTF = wlanLLTFDemodulate(lltf,cfgVHT.ChannelBandwidth);
        nVarVHT = helperNoiseEstimate(demodLLTF,cfgVHT);

        % Extract VHT-LTF samples from the waveform, demodulate and perform
        % channel estimation
        vhtltf = rx(ind.VHTLTF(1):ind.VHTLTF(2),:);
        vhtltfDemod = wlanVHTLTFDemodulate(vhtltf,cfgVHT);
        chanEst = wlanVHTLTFChannelEstimate(vhtltfDemod,cfgVHT);

        % Recover the transmitted PSDU in VHT Data
        % Extract VHT Data samples from the waveform and recover the PSDU
        vhtdata = rx(ind.VHTData(1):ind.VHTData(2),:);
        rxPSDU = wlanVHTDataRecover(vhtdata,chanEst,nVarVHT,cfgVHT);
        
        % Determine if any bits are in error, i.e. a packet error
        packetError = any(biterr(txPSDU,rxPSDU));
        numPacketErrors = numPacketErrors+packetError;
        numPkt = numPkt+1;
    end

    % Calculate packet error rate (PER) at SNR point
    packetErrorRate(i) = numPacketErrors/(numPkt-1);
    disp(['SNR ' num2str(snr(i)) ' completed after ' ...
        num2str(numPkt-1) ' packets, PER: ' ... 
        num2str(packetErrorRate(i))]);
end

%% Plot Packet Error Rate vs SNR Results

figure
semilogy(snr,packetErrorRate,'-ob');
grid on;
xlabel('SNR (dB)');
ylabel('PER');
title('802.11ac 80MHz, MCS9, Direct Mapping, 8x8 Channel Model D-NLOS');

%% Further Exploration
% The number of packets tested at each SNR point is controlled by two
% parameters; |maxNumErrors| and |maxNumPackets|. For meaningful results it
% is recommend that these values should be larger than those presented in
% this example. Increasing the number of packets simulated allows the PER
% under different scenarios to be compared. Try changing the transmission
% and reception configurations and compare the packet error rate. As an
% example, the figure below was created by running the example for longer.
%
% <<VHTMIMOPERExample.png>>

%% Appendix
% This example uses the following helper functions:
%
% * <matlab:edit('helperSampleRate.m') helperSampleRate.m>
% * <matlab:edit('helperFFTLength.m') helperFFTLength.m>
% * <matlab:edit('helperSubcarrierIndices.m') helperSubcarrierIndices.m>
% * <matlab:edit('helperPacketDetect.m') helperPacketDetect.m>
% * <matlab:edit('helperSymbolTiming.m') helperSymbolTiming.m>
% * <matlab:edit('helperNoiseEstimate.m') helperNoiseEstimate.m>

%% Selected Bibliography
% # IEEE Std 802.11ac(TM)-2013 IEEE Standard for Information technology -
% Telecommunications and information exchange between systems - Local and
% metropolitan area networks - Specific requirements - Part 11: Wireless
% LAN Medium Access Control (MAC) and Physical Layer (PHY) Specifications -
% Amendment 4: Enhancements for Very High Throughput for Operation in Bands
% below 6 GHz.

displayEndOfDemoMessage(mfilename)