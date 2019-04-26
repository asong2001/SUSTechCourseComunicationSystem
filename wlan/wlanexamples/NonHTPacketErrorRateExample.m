%% 802.11p(TM) and 802.11a(TM) Packet Error Rate Simulations 
%
% This example shows how to measure packet error rates of IEEE(R)
% 802.11p(TM) and 802.11a(TM) links using an end-to-end simulation with a
% fading channel model and additive white Gaussian noise. For similar link
% parameters, the two links are simulated alongside to offer a comparative
% perspective using the packet error rate as the figure of merit. The
% comparison highlights the robustness of the 802.11p link when compared
% with the 802.11a link.

% Copyright 2015 The MathWorks, Inc.

%% Introduction
% IEEE 802.11p [ <#10 2> ] is an approved amendment to the IEEE 802.11
% standard to enable support for wireless access in vehicular environments
% (WAVE). Using the half-clocked mode with a 10 MHz channel bandwidth, it
% operates at the 5.85-5.925 GHz bands to support applications for
% Intelligent Transportation Systems (ITS) [ <#10 3> ].
%
% In this example an end-to-end simulation is used to determine the packet
% error rates for 802.11p [ <#10 2> ] and 802.11a [ <#10 1> ] links with a
% fading channel at a selection of SNR points. For each SNR point multiple
% packets are transmitted through a fading channel, demodulated and the
% PSDUs recovered. The PSDUs are compared to those transmitted to determine
% the number of packet errors. Front-end components including packet
% detection, timing synchronization, carrier frequency offset correction
% and phase tracking are optionally enabled for the receiver. A schematic
% for the per-link processing is shown below.
%
% <<nonHTPERSchematic.png>>
%

%% Waveform Configuration
% Non-HT format transmissions are simulated in this example. For 802.11p, a
% 10 MHz channel bandwidth is used while for 802.11a a 20 MHz channel
% bandwidth is used. The two individual format parameters are specified
% using Non-HT format configuration objects by the
% <matlab:doc('wlanNonHTConfig') wlanNonHTConfig> function. In this example
% both links are configured for 16QAM rate 1/2 (MCS 4) operation.

% Link parameters
mcs = 4;                % QPSK rate 1/2
psduLen = 1000;         % PSDU length in bytes

% Create a format configuration object for a 802.11p transmission
cfgNHT10 = wlanNonHTConfig;
cfgNHT10.ChannelBandwidth = 'CBW10';    % 10 MHz channel bandwidth
cfgNHT10.PSDULength = psduLen;          
cfgNHT10.MCS = mcs;                     

% Create a format configuration object for a 802.11a transmission
cfgNHT20 = wlanNonHTConfig;
cfgNHT20.ChannelBandwidth = 'CBW20';    % 20 MHz channel bandwidth
cfgNHT20.PSDULength = psduLen;          
cfgNHT20.MCS = mcs;                     

%% Channel Configuration
% Both 802.11a and 802.11p are single antenna transmit systems. The
% difference in the channel bandwidth between the two links leads to a
% doubling of the signaling duration as seen by the following table which
% highlights the parametric differences for the two.
%
% <<nonHTPERparams.png>>
%
% The example uses a HIPERLAN/2 SISO fading channel model with delay
% profile Model-E representing NLOS conditions with an average delay spread
% of 250 ns [ <#10 3> ]. Such a channel with larger delays is applicable
% for outdoor vehicular operation. HIPERLAN/2 channel models with different
% delay profiles can be modeled using the <matlab:doc('stdchan') stdchan>
% function. For more representative 802.11a indoor scenarios, refer to
% the '802.11a' channel type option for the stdchan function.

% Create and configure the channel
chanMdl = 'E';
fd = 50;                            % Maximum Doppler shift, Hz
c = 3e8*3.6;                        % Speed of light, Km/hr
fc = 5.9e9;                         % Carrier frequency, Hz
disp(['Speed of unit = ' num2str(c*fd/5.9e9) ' Km/hr at ' num2str(fc/1e9) ' GHz']);

fs20 = helperSampleRate(cfgNHT20);     % Baseband sampling rate for 20 MHz
chan20 = stdchan(1/fs20, fd, ['hiperlan2' chanMdl]);

fs10 = helperSampleRate(cfgNHT10);     % Baseband sampling rate for 10 MHz
chan10 = stdchan(1/fs10, fd, ['hiperlan2' chanMdl]);

%% Simulation Parameters
% The operating SNR value per MCS value is simulated over a range of SNR
% points. For each SNR point simulated, a number of packets are generated,
% passed through the channel and demodulated to determine the packet error
% rate.

fadingMargin = 12; % dB
snrOperatingVec = [-1 1.75 2 4.75 7.5 10.75 15 16.5] + fadingMargin;
simRange = -2:1:1;
snr = snrOperatingVec(mcs+1) + simRange;

enableFE = false;    % Disable front-end receiver components

%% Simulation Setup
%
% Set up the simulation length per SNR point using a maximum number of
% packets and a maximum number of errors to be simulated, whichever occurs
% first.

maxNumErrors = 20;   % The maximum number of packet errors at an SNR point
maxNumPackets = 200; % Maximum number of packets at an SNR point

% Set random stream for repeatability of results
s = rng(98765);

%% Processing SNR Points
% For each SNR point a number of packets are tested and the packet error
% rate calculated for both links. For each packet the following processing
% steps occur:
%
% # A PSDU is created and encoded to create a single packet waveform.
% # The waveform is passed through a different realization of the 
% channel model.
% # AWGN is added to the received waveform to create the desired average
% SNR per subcarrier after OFDM demodulation.
% # Using the optional switch |enableFE|, front-end receiver components may
% be enabled. When enabled, the per-packet processing includes packet
% detection, coarse carrier frequency offset estimation and correction,
% symbol timing and fine carrier frequency offset estimation and
% correction. When disabled, the received waveform is synchronized using
% the known channel delay. 
% # The L-LTF is extracted from the synchronized received waveform. The
% L-LTF is OFDM demodulated and channel estimation using the L-LTF is
% performed.
% # The Non-HT Data field is extracted from the synchronized received
% waveform. The PSDU is recovered using the extracted data field and the
% L-LTF based channel estimate and noise power estimate.
%
% Refer to accompanying <matlab:edit('nonHTPERSimulator.m')
% nonHTPERSimulator.m> function for the processing details.

% Set up a figure for visualizing PER results
h = figure;
grid on;
hold on;
ax = gca;
ax.YScale = 'log';
xlim([snr(1), snr(end)]);
ylim([1e-3 1]);
xlabel('SNR (dB)');
ylabel('PER');
h.NumberTitle = 'off';
h.Name = '802.11p vs. 802.11a PER';
title(['MCS ' num2str(mcs) ', HIPERLAN/2 Model ' chanMdl ...
       ', Doppler ' num2str(fd) ' Hz']);

% Simulation loop for both links
S = numel(snr);
per20 = zeros(S,1); 
per10 = per20;
for i = 1:S     
    % 802.11p link
    per10(i) = nonHTPERSimulator(cfgNHT10, chan10, snr(i), ...
        maxNumErrors, maxNumPackets, enableFE);

    % 802.11a link
    per20(i) = nonHTPERSimulator(cfgNHT20, chan20, snr(i), ...
        maxNumErrors, maxNumPackets, enableFE);

    % Compare
    semilogy(snr, per10, 'bx-');
    semilogy(snr, per20, 'ro-');
    legend('802.11p, 10 MHz', '802.11a, 20 MHz');
    drawnow;
end
axis([5 35 1e-2 1])
hold off;

% Restore default stream
rng(s);

%% Conclusion 
%
% Observe the improved performance of the 10 MHz link for the higher delay
% spread channel, when compared with the similarly configured 20 MHz link.
% The longer cyclic prefix duration for 10 MHz operation mitigates the
% larger delay spread of the channel and offers increased robustness for
% outdoor operation.
%
% The simulation results presented are for a short run. For meaningful
% results, longer runs for larger number of packets and more errors are
% recommended. The result below is one such example.
%
% <<nonHTPERresults.png>>

%% Further Exploration
%
% This example highlights PER performance for only one of eight allowed MCS
% values. Explore link performance for other MCS values, varied Doppler
% shifts, different channel delay profiles and different receiver
% configurations.
%
% Alternate channel estimation and phase tracking algorithms that allow
% for higher speeds as needed by vehicular applications should also be
% considered.

%% Appendix
% This example uses the following helper functions:
%
% * <matlab:edit('helperSampleRate.m') helperSampleRate.m>
% * <matlab:edit('nonHTPERSimulator.m') nonHTPERSimulator.m>

%% Selected Bibliography
% # IEEE Std 802.11-2012: IEEE Standard for Information technology -
% Telecommunications and information exchange between systems - Local and
% metropolitan area networks - Specific requirements, Part 11: Wireless LAN
% Medium Access Control (MAC) and Physical Layer (PHY) Specifications,
% IEEE, New York, NY, USA, 1999-2013.
% # IEEE Std 802.11p-2010: IEEE Standard for Information technology -
% Telecommunications and information exchange between systems - Local and
% metropolitan area networks - Specific requirements, Part 11: Wireless LAN
% Medium Access Control (MAC) and Physical Layer (PHY) Specifications,
% Amendment 6: Wireless Access in Vehicular Environments,
% IEEE, New York, NY, USA, 2010.
% # ETSI,
% http://www.etsi.org/technologies-clusters/technologies/intelligent-transport/,
% 2015.
% # Medbo, J., P. Schramm, "Channel models for HIPERLAN/2", ETSI/BRAN 
% document no. 3ERI085B, 1998. 

displayEndOfDemoMessage(mfilename)