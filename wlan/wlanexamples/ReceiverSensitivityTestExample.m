%% 802.11ac(TM) Receiver Minimum Input Sensitivity Test
% This example shows how to simulate a test to measure the receiver minimum
% input sensitivity as specified in Section 22.3.19.1 of the IEEE(R)
% 802.11ac standard [ <#17 1> ].

% Copyright 2015 The MathWorks, Inc.

%% Introduction
% The receiver minimum sensitivity test ensures a device under test (DUT)
% is able to receive data with a defined maximum packet error rate (PER) of
% 10% at a defined minimum signal power. The minimum signal power depends
% on the channel bandwidth and modulation and coding scheme (MCS):
%
% <<VHTRxSensitivityMinimumTable.png>>
%
% When the test is performed with hardware, each input antenna port on the
% DUT is connected through a cable to a single output antenna port of a
% transmitter. The following transmission parameters are specified for the
% test waveform:
%
% * The number of spatial streams is equal to the number of transmitting
% antenna ports
% * PSDU length of 4096 bytes
% * No space-time block coding (STBC)
% * 800ns guard interval
% * Binary convolutional coding
%
% This example shows how the above test can be constructed with an
% end-to-end simulation using WLAN System Toolbox(TM). In this example a
% receiver is stimulated with incoming VHT packets at a range of input
% levels below the minimum sensitivity level and the PER measured.
%
% For each sensitivity level tested packets are generated and scaled to the
% desired signal level. Additive white Gaussian noise (AWGN) is added to
% create a noise floor at the receiver. The noisy packets are then
% demodulated and the PSDUs recovered. The PSDUs are compared to those
% transmitted to determine the number of packet errors and hence the packet
% error rate. Packet detection, timing synchronization, carrier frequency
% offset correction, noise estimation and phase tracking are performed by
% the receiver. The processing for each packet is summarized in the
% following diagram:
%
% <<VHTRxSensitivityDiagram.png>>

%% Test Parameters
% The transmission configuration for the test is specified with a VHT
% configuration object. In this example the minimum sensitivity is measured
% for a 160 MHz transmission with 64-QAM rate 5/6 modulation and coding.
% The simulated DUT has 2 receive antennas. These parameters can be changed
% to test different configurations.

cfgVHT = wlanVHTConfig;             % Create VHT transmission configuration
cfgVHT.ChannelBandwidth = 'CBW160'; % Bandwidth
cfgVHT.MCS = 7;                     % 64-QAM, rate 5/6
NumReceiveAntennas = 2;             % Number of receive antennas

%%
% The fixed transmission parameters required by the test are set below.

cfgVHT.APEPLength = 4096; % Bytes
cfgVHT.STBC = false;
cfgVHT.NumTransmitAntennas = NumReceiveAntennas;
cfgVHT.NumSpaceTimeStreams = NumReceiveAntennas;
cfgVHT.SpatialMapping = 'Direct';
cfgVHT.GuardInterval = 'Long';

%% Simulation Parameters
% In this example a receiver is stimulated with VHT packets at a range of
% input levels below the minimum input sensitivity level. The range of
% offsets tested is specified in the vector |testInputLevelOffsets|.

testInputLevelOffsets = [-10 -9 -8 -7]; % dB

%%
% The number of packets tested at each sensitivity is controlled by two
% parameters:
%
% # |maxNumErrors| is the maximum number of packet errors simulated at each
% input level. When the number of packet errors reaches this limit, the
% simulation at this sensitivity is complete.
% # |maxNumPackets| is the maximum number of packets simulated at each
% input level point and limits the length of the simulation if the packet
% error limit is not reached.
%
% The numbers chosen in this example will lead to a very short simulation.
% For meaningful results we recommend increasing the numbers.

maxNumErrors = 20;
maxNumPackets = 200;

%% Signal Power Setup
% The minimum sensitivity test specifies a maximum PER for a measured input
% level per receive antenna. In this simulation the receiver is stimulated
% with a test signal with a specified input level in dBm. The test signal
% is generated with the waveform generator,
% <matlab:doc('wlanWaveformGenerator') wlanWaveformGenerator>. The output
% of the waveform generator is normalized internally such that the sum of
% the power for all antenna powers is 0 dBm. Therefore for this simulation
% the output of the waveform generator must be scaled to create the desired
% input level.
%
% First the minimum sensitivity for the transmission configuration is
% determined from Table 22-25 of the 802.11ac standard [ <#17 1> ].

% Receiver minimum input level sensitivity for 20 MHz, Table 22-25 Std
% 802.11-2013ac. The sensitivity increases by 3dB for double the bandwidth.
rxMinSensitivityTable = [-82 -79 -77 -74 -70 -66 -65 -64 -59 -57]; % dBm

% Get minimum input sensitivity given MCS and bandwidth
fs = helperSampleRate(cfgVHT);     % Baseband sampling rate (Hz)
B = floor(10*log10((fs/20e6))); % Scalar for bandwidth
rxMinSensitivity = rxMinSensitivityTable(cfgVHT.MCS+1)+B; % dBm
disp(['Minimum sensitivity for MCS' num2str(cfgVHT.MCS) ', ' ...
    num2str(fs/1e6) ' MHz: ' num2str(rxMinSensitivity,'%2.1f') ' dBm'])

%%
% In this example a range of input levels below the minimum level are
% tested. These power levels are defined by |testInputLevels|.

testInputLevels = rxMinSensitivity+testInputLevelOffsets; % dBm

%%
% For a given required test input level per antenna, a voltage scalar, |A|,
% to apply to the generated waveform is calculated. The power per receive
% antenna port is measured during the simulation to confirm the input
% signal level is correct.

A = 10.^((testInputLevels-30)/20);      % Voltage gain (attenuation)
A = A*sqrt(cfgVHT.NumTransmitAntennas); % Account for generator scaling

%% Noise Configuration
% The noise floor of the receiver is simulated with thermal noise. The
% height of the noise floor determines the SNR at the receiver, as the
% input signal level is fixed for this test. The noise figure of the
% receiver determines the level of noise floor.

NF = 6;         % Noise figure (dB)
T = 290;        % Ambient temperature (K)
BW = fs;        % Bandwidth (Hz)
k = 1.3806e-23; % Boltzmann constant (J/K)
noiseFloor = 10*log10(k*T*BW)+NF; % dB
disp(['Receiver noise floor: ' num2str(noiseFloor+30,'%2.1f') ' dBm'])

%%
% An AWGN channel, <matlab:doc('comm.AWGNChannel') comm.AWGNChannel>, is
% used to add noise to the waveform.

AWGN = comm.AWGNChannel('NoiseMethod','Variance', ...
    'Variance',10^(noiseFloor/10));

%% Simulation Setup
% In this simulation frequency offset estimation is used.
% <matlab:doc('comm.PhaseFrequencyOffset') comm.PhaseFrequencyOffset>
%  is used to correct the carrier frequency offset estimated
% by <matlab:doc('wlanCoarseCFOEstimate') wlanCoarseCFOEstimate> and
% <matlab:doc('wlanFineCFOEstimate') wlanFineCFOEstimate>.

PFO = comm.PhaseFrequencyOffset('FrequencyOffsetSource','Input port', ...
    'SampleRate',fs);

%%
% Set the remaining parameters for the simulation.

% Indices for accessing each field within the time-domain packet
ind = wlanFieldIndices(cfgVHT);

chanBW = cfgVHT.ChannelBandwidth;

rng(0);  % Set random state for repeatability

%% Input Level Sensitivity Simulation
% For each input level a number of packets are tested and the packet error
% rate calculated.
%
% For each packet the following processing steps occur:
%
% # A PSDU is created and encoded to create a single packet waveform.
% # The waveform is scaled to create the desired input level in dBm.
% # The power of the received waveform is measured.
% # AWGN is added to the received waveform to create a noise floor.
% # The packet is detected.
% # Coarse carrier frequency offset is estimated and corrected.
% # Fine timing synchronization is established. The L-STF, L-LTF and L-SIG
% samples are provided for fine timing to allow for packet detection at the
% start or end of the L-STF.
% # Fine carrier frequency offset is estimated and corrected.
% # The noise power is estimated using the L-LTF.
% # The VHT-LTF is extracted from the synchronized received waveform. The
% VHT-LTF is OFDM demodulated and channel estimation is performed.
% # The VHT Data field is extracted from the synchronized received
% waveform. The PSDU is recovered using the extracted field and the channel
% estimate.

S = numel(testInputLevels);
packetErrorRate = zeros(S,1);
rxAntennaPower = zeros(S,1);
for i=1:S
    disp(['Simulating ' num2str(testInputLevels(i),'%2.1f') ...
        ' dBm input level...']);
    
    % Loop to simulate multiple packets
    numPacketErrors = 0;
    measuredPower = []; % Measured average power per antenna
    numPkt = 1; % Index of packet transmitted
    while numPacketErrors<=maxNumErrors && numPkt<=maxNumPackets
        % Generate a packet waveform
        txPSDU = randi([0 1],cfgVHT.PSDULength*8,1); % PSDULength in bytes
        tx = wlanWaveformGenerator(txPSDU,cfgVHT);
        
        % Scale input signal to desired level
        rx = tx.*A(i);

        % Measure the average power at the antenna connector in Watts
        measuredPower(numPkt) = mean(mean(rx.*conj(rx))); %#ok<SAGROW>

        % Add noise floor at receiver
        rx = step(AWGN,rx);
        
        % Packet detect
        pktStartIdx = helperPacketDetect(rx,chanBW);
        if isempty(pktStartIdx) % If empty no L-STF detected; packet error
            numPacketErrors = numPacketErrors+1;
            numPkt = numPkt+1;
            continue; % Go to next loop iteration
        end
        pktOffset = pktStartIdx-1; % Packet offset from start of waveform

        % Extract L-STF and perform coarse frequency offset correction
        lstf = rx(pktOffset+(ind.LSTF(1):ind.LSTF(2)),:); 
        coarseFreqOff = wlanCoarseCFOEstimate(lstf,chanBW);
        rx = step(PFO,rx,-coarseFreqOff);
        release(PFO); % Release object for subsequent processing

        % Extract the Non-HT fields and determine start of L-LTF
        nonhtfields = rx(pktOffset+(ind.LSTF(1):ind.LSIG(2)),:); 
        lltfIdx = helperSymbolTiming(nonhtfields,chanBW);

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
        fineFreqOff = wlanFineCFOEstimate(lltf,chanBW);
        rx = step(PFO,rx,-fineFreqOff);
        release(PFO); % Release object for subsequent processing

        % Estimate noise power in VHT fields
        lltf = rx(ind.LLTF(1):ind.LLTF(2),:);
        demodLLTF = wlanLLTFDemodulate(lltf,cfgVHT);
        nEstVHT = helperNoiseEstimate(demodLLTF,cfgVHT);
                
        % Extract VHT-LTF samples from the waveform, demodulate and perform
        % channel estimation
        vhtltf = rx(ind.VHTLTF(1):ind.VHTLTF(2),:);
        vhtltfDemod = wlanVHTLTFDemodulate(vhtltf,cfgVHT);
        chanEst = wlanVHTLTFChannelEstimate(vhtltfDemod,cfgVHT);

        % Recover the transmitted PSDU in VHT Data
        % Extract VHT Data samples from the waveform and recover the PSDU
        vhtdata = rx(ind.VHTData(1):ind.VHTData(2),:);
        rxPSDU = wlanVHTDataRecover(vhtdata,chanEst,nEstVHT,cfgVHT);
        
        % Determine if any bits are in error, i.e. a packet error
        packetError = any(biterr(txPSDU,rxPSDU));
        numPacketErrors = numPacketErrors+packetError;
        numPkt = numPkt+1;
    end

    % Calculate packet error rate (PER) at input level point
    packetErrorRate(i) = numPacketErrors/(numPkt-1);
    disp(['  Completed after ' ...
        num2str(numPkt-1) ' packets, PER: ' ... 
        num2str(packetErrorRate(i))]);

    % Calculate average input power per antenna
    rxAntennaPower(i) = 10*log10(mean(measuredPower))+30;
    disp(['  Measured antenna connector power: ' ...
        num2str(rxAntennaPower(i),'%2.1f') ' dBm']);
end

%% Analysis and Further Exploration
% The PER for tested input signal levels is plotted with the maximum PER at
% minimum sensitivity.

figure
semilogy(rxAntennaPower,packetErrorRate,'o-')
hold on
semilogy(rxMinSensitivity,0.1,'rx')
currentxlim = xlim(gca);
xlim([currentxlim(1) currentxlim(2)+1])
grid on
xlabel('Measured power per antenna connector (dBm)');
ylabel('PER');
legend('Simulated PER performance','Maximum PER at minimum sensitivity');
title(sprintf(['Minimum Input Sensitivity Test: MCS%d, %d MHz, ' ...
    '%d Antennas'],cfgVHT.MCS,fs/1e6,cfgVHT.NumTransmitAntennas))

%%
% Inspecting the plot reveals the simulated 10% PER is just under 8 dB
% lower than the minimum sensitivity specified by the test. This difference
% is due to the implementation margin allowed by the test. The
% implementation margin allows for algorithmic degradations due to
% impairments and the receiver noise figure when compared to ideal AWGN
% performance [ <#17 2> ]. In this example only AWGN is added as an
% impairment. Therefore only the algorithmic performance of front-end
% synchronization, channel estimation and phase tracking in the presence of
% AWGN use the implementation margin. If more impairments are included in
% the simulation the PER waterfall in the plot will move right towards the
% minimum sensitivity and the margin will decrease.
%
% The number of packets tested at each SNR point is controlled by two
% parameters; |maxNumErrors| and |maxNumPackets|. For meaningful results it
% is recommend that these values should be larger than those presented in
% this example. 

%% Appendix
% This example uses the following helper functions:
%
% * <matlab:edit('helperSampleRate.m') helperSampleRate.m>
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
% # Perahia, Eldad, and Robert Stacey. Next Generation Wireless LANS:
% 802.11n and 802.11ac. Cambridge University Press, 2013.

displayEndOfDemoMessage(mfilename)