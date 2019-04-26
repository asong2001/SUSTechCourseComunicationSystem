%% 802.11n(TM) Packet Error Rate Simulation for 2x2 TGn Channel
%
% This example shows how to measure the packet error rate of an IEEE(R)
% 802.11n(TM) HT link using an end-to-end simulation with a fading TGn
% channel model and additive white Gaussian noise.

% Copyright 2015 The MathWorks, Inc.

%% Introduction
% In this example an end-to-end simulation is used to determine the packet
% error rate for an 802.11n HT [ <#12 1> ] link with a fading channel at a
% selection of SNR points. At each SNR point multiple packets are
% transmitted through a channel, demodulated and the PSDUs recovered. The
% PSDUs are compared to those transmitted to determine the number of packet
% errors and hence the packet error rate. Packet detection, timing
% synchronization, carrier frequency offset correction and phase tracking
% are performed by the receiver. The processing for each packet is
% summarized in the following diagram.
%
% <<HTMIMOPERDiagram.png>>
%
% This example also demonstrates how a <matlab:doc('parfor') parfor> loop
% can be used instead of the <matlab:doc('for') for> loop when simulating
% each SNR point to speed up a simulation. <matlab:doc('parfor') parfor>,
% as part of the Parallel Computing Toolbox(TM), executes processing for
% each SNR in parallel to reduce the total simulation time.

%% Waveform Configuration
% An 802.11n HT transmission is simulated in this example. The HT format
% configuration object contains the format specific configuration of the
% transmission. The object is created using the <matlab:doc('wlanHTConfig')
% wlanHTConfig> function. The properties of the object contain the
% configuration. In this example the object is configured for a 20 MHz
% channel bandwidth, 2 transmit antennas, 2 space time streams and no space
% time block coding.

% Create a format configuration object for a 2-by-2 HT transmission
cfgHT = wlanHTConfig;
cfgHT.ChannelBandwidth = 'CBW20'; % 20 MHz channel bandwidth
cfgHT.NumTransmitAntennas = 2;    % 2 transmit antennas
cfgHT.NumSpaceTimeStreams = 2;    % 2 space-time streams
cfgHT.PSDULength = 1000;          % PSDU length in bytes
cfgHT.MCS = 15;                   % 2 spatial streams, 64-QAM rate-5/6

%% Channel Configuration
% In this example a TGn N-LOS channel model is used with delay profile
% Model-B. For Model-B when the distance between transmitter and receiver
% is greater than or equal to five meters, the model is NLOS. This is
% described further in <matlab:doc('wlanTGnChannel') wlanTGnChannel>.

% Create and configure the channel
tgn = wlanTGnChannel;
tgn.DelayProfile = 'Model-B';
tgn.NumTransmitAntennas = cfgHT.NumTransmitAntennas;
tgn.NumReceiveAntennas = 2;
tgn.TransmitReceiveDistance = 10; % Distance in meters for NLOS operation
tgn.LargeScaleFadingEffect = 'None';

%% Simulation Parameters
% For each SNR point in the vector |snr| a number of packets are
% generated, passed through a channel and demodulated to determine the
% packet error rate.

snr = 25:10:45;

%%
% The number of packets tested at each SNR point is controlled by two
% parameters:
%
% # |maxNumPEs| is the maximum number of packet errors simulated at each
% SNR point. When the number of packet errors reaches this limit, the
% simulation at this SNR point is complete.
% # |maxNumPackets| is the maximum number of packets simulated at each SNR
% point and limits the length of the simulation if the packet error limit
% is not reached. 
%
% The numbers chosen in this example will lead to a very short simulation.
% For meaningful results we recommend increasing the numbers.

maxNumPEs = 10; % The maximum number of packet errors at an SNR point
maxNumPackets = 100; % Maximum number of packets at an SNR point

%% Simulation Setup
% In this simulation frequency offset estimation is used.
% <matlab:doc('comm.PhaseFrequencyOffset') comm.PhaseFrequencyOffset>
%  is used to correct the carrier frequency offset estimated
% by <matlab:doc('wlanCoarseCFOEstimate') wlanCoarseCFOEstimate> and
% <matlab:doc('wlanFineCFOEstimate') wlanFineCFOEstimate>.

% Get the baseband sampling rate
fs = helperSampleRate(cfgHT);

% Create an instance of the frequency impairment object
PFO = comm.PhaseFrequencyOffset;
PFO.SampleRate = fs;
PFO.PhaseOffset = 0;
PFO.FrequencyOffsetSource = 'Input port';

%%
% Set the remaining parameters for the simulation.

% Get the number of occupied subcarriers in HT fields and FFT length
[htData,htPilots] = helperSubcarrierIndices(cfgHT,'HT');
Nst_ht = numel(htData)+numel(htPilots);
Nfft = helperFFTLength(cfgHT);      % FFT length

% Set the sampling rate of the channel
tgn.SampleRate = fs;

% Indices for accessing each field within the time-domain packet
ind = wlanFieldIndices(cfgHT);

rng(0); % Set random state for repeatability

%% Processing SNR Points
% For each SNR point a number of packets are tested and the packet error
% rate calculated.
%
% For each packet the following processing steps occur:
%
% # A PSDU is created and encoded to create a single packet waveform.
% # The waveform is passed through a different realization of the TGn
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
% # The HT-LTF is extracted from the synchronized received waveform. The
% HT-LTF is OFDM demodulated and channel estimation is performed.
% # The HT Data field is extracted from the synchronized received waveform.
% The PSDU is recovered using the extracted field and the channel estimate.
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
for i = 1:S % Use 'for' to debug the simulation
    % Create an instance of the AWGN channel per SNR point simulated
    AWGN = comm.AWGNChannel;
    AWGN.NoiseMethod = 'Signal to noise ratio (SNR)';
    AWGN.SignalPower = 1/tgn.NumReceiveAntennas; % Normalization
    AWGN.SNR = snr(i)-10*log10(Nfft/Nst_ht); % Account for energy in nulls

    % Loop to simulate multiple packets
    numPacketErrors = 0;
    n = 1; % Index of packet transmitted
    while numPacketErrors<=maxNumPEs && n<=maxNumPackets
        % Generate a packet waveform
        txPSDU = randi([0 1],cfgHT.PSDULength*8,1); % PSDULength in bytes
        tx = wlanWaveformGenerator(txPSDU,cfgHT);
        
        % Add trailing zeros to allow for channel filter delay
        tx = [tx; zeros(15,cfgHT.NumTransmitAntennas)]; %#ok<AGROW>
        
        % Pass the waveform through the TGn channel model 
        rx = step(tgn,tx);
        reset(tgn); % Reset channel to create a different realization 
                
        % Add noise
        rx = step(AWGN,rx);
 
        % Packet detect
        pktStartIdx = helperPacketDetect(rx,cfgHT.ChannelBandwidth);
        if isempty(pktStartIdx) % If empty no L-STF detected; packet error
            numPacketErrors = numPacketErrors+1;
            n = n+1;
            continue; % Go to next loop iteration
        end
        pktOffset = pktStartIdx-1; % Packet offset from start of waveform
        
        % Extract L-STF and perform coarse frequency offset correction
        lstf = rx(pktOffset+(ind.LSTF(1):ind.LSTF(2)),:); 
        coarseFreqOff = wlanCoarseCFOEstimate(lstf,cfgHT.ChannelBandwidth);
        rx = step(PFO,rx,-coarseFreqOff);
        release(PFO); % Release object for subsequent processing
        
        % Extract the Non-HT fields and determine start of L-LTF
        nonhtfields = rx(pktOffset+(ind.LSTF(1):ind.LSIG(2)),:); 
        lltfIdx = helperSymbolTiming(nonhtfields,cfgHT.ChannelBandwidth);
        
        % Synchronize the received waveform given the offset between the
        % expected start of the L-LTF and actual start of L-LTF
        pktOffset = pktOffset+lltfIdx-double(ind.LLTF(1));
        % If no L-LTF detected or if packet detected outwith the range of
        % expected delays from the channel modeling; packet error
        if isempty(lltfIdx) || pktOffset<0 || pktOffset>15
            numPacketErrors = numPacketErrors+1;
            n = n+1;
            continue; % Go to next loop iteration
        end  
        rx = rx(1+pktOffset:end,:);
             
        % Extract L-LTF and perform fine frequency offset correction
        lltf = rx(ind.LLTF(1):ind.LLTF(2),:); 
        fineFreqOff = wlanFineCFOEstimate(lltf,cfgHT.ChannelBandwidth);
        rx = step(PFO,rx,-fineFreqOff);
        release(PFO); % Release object for subsequent processing
        
        % Estimate noise power in HT fields
        lltf = rx(ind.LLTF(1):ind.LLTF(2),:);
        demodLLTF = wlanLLTFDemodulate(lltf,cfgHT.ChannelBandwidth);
        nVarHT = helperNoiseEstimate(demodLLTF,cfgHT);
        
        % Extract HT-LTF samples from the waveform, demodulate and perform
        % channel estimation
        htltf = rx(ind.HTLTF(1):ind.HTLTF(2),:);
        htltfDemod = wlanHTLTFDemodulate(htltf,cfgHT);
        chanEst = wlanHTLTFChannelEstimate(htltfDemod,cfgHT);

        % Recover the transmitted PSDU in HT Data
        % Extract HT Data samples from the waveform and recover the PSDU
        htdata = rx(ind.HTData(1):ind.HTData(2),:);
        rxPSDU = wlanHTDataRecover(htdata,chanEst,nVarHT,cfgHT);
        
        % Determine if any bits are in error, i.e. a packet error
        packetError = any(biterr(txPSDU,rxPSDU));
        numPacketErrors = numPacketErrors+packetError;
        n = n+1;
    end
    
    % Calculate packet error rate (PER) at SNR point
    packetErrorRate(i) = numPacketErrors/(n-1);
    disp(['SNR ' num2str(snr(i)) ' completed after ' ...
        num2str(n-1) ' packets']);
end

%% Plot Packet Error Rate vs SNR Results
figure;
semilogy(snr,packetErrorRate,'-ob');
grid on;
xlabel('SNR [dB]');
ylabel('PER');
title('802.11n 20MHz, MCS15, Direct Mapping, 2x2 Channel Model B-NLOS');

%% Further Exploration
% The number of packets tested at each SNR point is controlled by two
% parameters; |maxNumPEs| and |maxNumPackets|. For meaningful results it is
% recommended that these values should be larger than those presented in
% this example. Increasing the number of packets simulated allows the PER
% under different scenarios to be compared. Try changing the transmission
% and reception configurations and compare the packet error rate. As an
% example, the figure below was created by running the example for four
% different configurations; 1x1, 2x2, 3x3 and 4x4.
%
% <<HTMIMOPERExample.png>>

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
% # IEEE Std 802.11(TM)-2012 IEEE Standard for Information technology -
% Telecommunications and information exchange between systems - Local and
% metropolitan area networks - Specific requirements - Part 11: Wireless
% LAN Medium Access Control (MAC) and Physical Layer (PHY) Specifications.

displayEndOfDemoMessage(mfilename)