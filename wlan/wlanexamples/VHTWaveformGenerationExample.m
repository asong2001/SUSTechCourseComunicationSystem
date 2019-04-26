%% 802.11ac(TM) Waveform Generation with MAC Frames
%
% This example shows how to generate an IEEE(R) 802.11ac(TM) transmission
% containing MAC frames suitable for performing radio packet error rate
% (PER) receiver tests.

% Copyright 2015 The MathWorks, Inc.

%% Introduction
% WLAN System Toolbox(TM) can be used to generate standard compliant
% waveforms for performing receiver tests. A basic WLAN receiver test
% scenario is shown in the diagram below.
%
% <<vhtWaveformGenerationDiagram.png>>
%
% The device under test (DUT) is stimulated with RF test vectors, usually
% through a wired link. The packet error rate (PER) is a metric used to
% test the performance of a receiver at a given receive signal power in
% the presence of noise, interference, or other impairments. The PER is
% defined as the number of incorrectly decoded packets divided by the total
% number of transmitted packets.
%
% The frame check sum (FCS) within a MAC frame is used to determine whether
% a MAC frame has been decoded correctly by the receiver, and therefore
% whether the packet has been received in error. The general MAC frame for
% IEEE 802.11ac contains the following fields:
% 
% * MAC header
% * Frame body
% * FCS
%
% The data to transmit from a higher layer is contained within the frame
% body of the MAC frame. The transmitter uses a cyclic redundancy check
% over the MAC header and frame body field to generate the FCS value. The
% receiver calculates the CRC and compares this to the received FCS field
% to determine if an error has occurred during transmission.
% 
% In this example an IEEE 802.11ac waveform consisting of multiple VHT
% format packets is generated. The function
% <matlab:doc('wlanWaveformGenerator') wlanWaveformGenerator> can be used
% to generate a waveform containing one or more packets.
% <matlab:doc('wlanWaveformGenerator') wlanWaveformGenerator> consumes
% physical layer service data units (PSDUs) for each packet and performs
% the appropriate physical layer processing to create the waveform. In this
% example a multi-packet baseband waveform with a MAC header and valid FCS
% is synthesized. This waveform may be downloaded to a signal generator for
% RF transmission and used for receiver PER testing. Source code is
% provided to download and play the waveform using a Keysight Technologies(TM) 
% N5172B signal generator. The example processing is illustrated in the
% following diagram:
%
% <<vhtWaveformExampleDiagram.png>>

%% IEEE 802.11ac VHT Format Configuration
% The format specific configuration of a VHT waveform synthesized with
% <matlab:doc('wlanWaveformGenerator') wlanWaveformGenerator> is described
% by a VHT format configuration object. The object is created using the
% <matlab:doc('wlanVHTConfig') wlanVHTConfig> function. The properties of
% the object contain the configuration. In this example an object is
% configured for a 160 MHz bandwidth, 1 transmit antenna, 1 space-time
% stream and QPSK rate 1/2 (MCS 1).

vhtCfg = wlanVHTConfig;             % Create packet configuration
vhtCfg.ChannelBandwidth = 'CBW160'; % 160 MHz channel bandwidth
vhtCfg.NumTransmitAntennas = 1;     % 1 transmit antenna
vhtCfg.NumSpaceTimeStreams = 1;     % 1 space-time stream
vhtCfg.MCS = 1;                     % Modulation: QPSK Rate: 1/2

%% Waveform Generation Configuration 
% The waveform generator <matlab:doc('wlanWaveformGenerator')
% wlanWaveformGenerator> can be configured to generate one or more packets
% and add an idle time between each packet. In this example four packets
% with a 20 microsecond idle period will be created.

numPackets = 4;   % Generate 4 packets
idleTime = 20e-6; % 20 microseconds idle period after packet

%%
% The PSDU transmitted in each packet is scrambled using a random seed for
% each packet. This is accomplished by specifying a vector of scrambler
% initialization seeds. The valid range of the seed is between 1 and 127
% inclusive.

% Initialize the scrambler with a random integer for each packet
scramblerInitialization = randi([1 127],numPackets,1);

%% Create a PSDU Containing a MAC Header, Body and FCS for each Packet
% For an IEEE 802.11ac data transmission the MAC frame is termed a MAC
% protocol data unit (MPDU), the MAC header is termed the MPDU header, and
% the frame body is an aggregated MAC service data unit (A-MSDU). One or
% more MPDUs are delimited, padded and aggregated to create an aggregated
% MPDU (A-MPDU). The A-MPDU is delimited and padded to form the physical
% layer service data unit (PSDU) which is coded and modulated to create the
% transmitted packet. This process of encapsulation is shown in the
% following diagram:
%
% <<vhtWaveformGenerationEncapsulation.png>>
%
% In this example a PSDU is created containing a single MPDU for each
% packet. The MPDU consists of an MPDU header containing no content, random
% frame body and valid FCS. The helper function
% <matlab:edit('vhtWaveformMACHeader.m') vhtWaveformMACHeader.m> creates a
% MPDU header with no content. A CRC is used to generate the FCS for the
% A-MPDU. The 32 bit long FCS is appended to the MPDU header and frame
% body.
%
% The helper function <matlab:edit('vhtWaveformGeneratePSDU.m')
% vhtWaveformGeneratePSDU.m> pads and delimits the MPDU to form a PSDU
% ready for physical layer processing and transmission by
% <matlab:doc('wlanWaveformGenerator') wlanWaveformGenerator> as specified
% in IEEE Std 802.11ac-2013 [ <#11 1> ]. The MPDU bits are delimited and
% form an A-MPDU containing a single MPDU. The length of this A-MPDU in
% bytes is the APEP Length, and the appropriate |vhtCfg.APEPLength|
% property of the VHT configuration object is set accordingly. The A-MPDU
% is then padded to form a PSDU of the required length, ready for physical
% layer processing. A PSDU is generated for each packet and is concatenated
% into a vector |data| for transmission with
% <matlab:doc('wlanWaveformGenerator') wlanWaveformGenerator>. The
% processing to create the concatenated PSDU bits |data| is shown in the
% diagram below:
%
% <<vhtWaveformMPDUProcessing.png>>

% Generate FCS and form MPDU, IEEE Std 802.11-2012 Section 8.2.4
fcsGen = comm.CRCGenerator([32 26 23 22 16 12 11 10 8 7 5 4 2 1 0]);
fcsGen.InitialConditions = 1;
fcsGen.DirectMethod = true;
fcsGen.FinalXOR = 1;

data = [];
for i=1:numPackets
    % Generate a payload for the frame body, a random A-MSDU of 4048 octets
    amsdu = randi([0 1],4048*8,1);
    mpduHeader = vhtWaveformMACHeader;

    % Create MPDU with header, body and FCS
    mpdu = step(fcsGen,[mpduHeader; amsdu]);
    
    % Delimit and pad the MPDU to form the PSDU for a single packet
    [psdu,vhtCfg] = vhtWaveformGeneratePSDU(mpdu,vhtCfg);
    
    % Concatenate packet PSDUs for waveform generation
    data = [data; psdu]; %#ok<AGROW>
end

%% Generate a Baseband Waveform
% The concatenated PSDU bits for all packets, |data|, are passed as an
% argument to <matlab:doc('wlanWaveformGenerator') wlanWaveformGenerator>
% along with the VHT packet configuration object |vhtCfg|. This configures
% the waveform generator to synthesize an 802.11ac VHT waveform. To
% generate 802.11n HT or other format waveforms different configuration
% objects are used such as those created by <matlab:doc('wlanHTConfig')
% wlanHTConfig> or <matlab:doc('wlanNonHTConfig') wlanNonHTConfig>. The
% waveform generator is additionally configured using name-value pairs to
% generate multiple packets with a specified idle time between packets, and
% initial scrambler states.

% Generate baseband VHT packets
txWaveform = wlanWaveformGenerator(data,vhtCfg, ...
    'NumPackets',numPackets,'IdleTime',idleTime, ...
    'ScramblerInitialization',scramblerInitialization);

fs = helperSampleRate(vhtCfg);
disp(['Baseband sampling rate: ' num2str(fs/1e6) ' Msps']);

%%
% The magnitude of the baseband waveform is displayed below. Note the
% number of packets and idle time configured.

figure;
plot(abs(txWaveform));
xlabel('Sample index');
ylabel('Magnitude');
title('Baseband IEEE 802.11ac Waveform');
legend('Transmit antenna 1');

%%
% The frequency spectrum of the generated time domain waveform,
% |txWaveform|, can be viewed using the
% <http://www.mathworks.com/products/dsp-system/ DSP System Toolbox(TM)>
% <matlab:doc('dsp.SpectrumAnalyzer') spectrum analyzer>. As expected, the
% 160 MHz signal bandwidth is clearly visible at baseband.

sa = dsp.SpectrumAnalyzer;
sa.SampleRate = fs;
sa.SpectrumType = 'Power density';
sa.RBWSource = 'Property';
sa.RBW = 100e3;
sa.FrequencySpan = 'Span and center frequency';
sa.Span = fs;
sa.SpectralAverages = 10;
sa.YLabel = 'PSD';
sa.YLimits = [-80 -40];
sa.Title = 'Baseband IEEE 802.11ac Waveform';
step(sa,txWaveform);

%% Generate an Over-the-Air Signal using an RF Signal Generator
% The baseband waveform created by WLAN System Toolbox can now be
% downloaded to a signal generator to perform receiver tests.
% <http://www.mathworks.com/products/instrument/ Instrument Control Toolbox(TM)> 
% is used to generate an RF signal with a center frequency of 5.25 GHz RF
% using the Keysight Technologies N5172B signal generator.

% Control whether to download the waveform to the waveform generator
playOverTheAir = false;

% Download the baseband IQ waveform to the instrument. Generate the RF 
% signal at a center frequency of 5.25 GHz and output power of -10 dBm.
if playOverTheAir
    fc = 5.25e9; %#ok<UNRCH> % Center frequency
    power = -10; % Output power
    hDownloadAndPlayWaveformUsingN5172B('192.168.0.1',txWaveform, ...
        fs,fc,power);
end

%% Appendix
% This example uses the following helper functions:
%
% * <matlab:edit('vhtWaveformMACHeader.m') vhtWaveformMACHeader.m>
% * <matlab:edit('vhtWaveformGeneratePSDU.m') vhtWaveformGeneratePSDU.m>
% * <matlab:edit('hDownloadAndPlayWaveformUsingN5172B.m') hDownloadAndPlayWaveformUsingN5172B.m>
% * <matlab:edit('helperSampleRate.m') helperSampleRate.m>

%% Selected Bibliography
% # IEEE Std 802.11ac(TM)-2013 IEEE Standard for Information technology -
% Telecommunications and information exchange between systems - Local and
% metropolitan area networks - Specific requirements - Part 11: Wireless
% LAN Medium Access Control (MAC) and Physical Layer (PHY) Specifications -
% Amendment 4: Enhancements for Very High Throughput for Operation in Bands
% below 6 GHz.

displayEndOfDemoMessage(mfilename)
