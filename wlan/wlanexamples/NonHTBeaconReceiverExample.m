%% 802.11(TM) OFDM Beacon Receiver with Live Data
%
% This example shows a receiver design that is able to recover 802.11(TM)
% OFDM Non-HT based beacon packets transmitted over the air from commercial
% 802.11 hardware. Beacon packets are typically transmitted in Non-HT
% format, even for HT [ <#9 1> ] and/or VHT [ <#9 2> ] capable hardware.
% Packet information such as SSID is printed to the command-line during
% recovery.

% Copyright 2015 The MathWorks, Inc.

%% Introduction
% This example illustrates the use of WLAN System Toolbox(TM) to recover
% real-world signals. It demonstrates a receiver design including
% synchronization, transmission configuration recovery, and payload
% decoding for Non-HT packets. Although this example recovers beacon
% packets from a file containing a captured baseband waveform, it is
% designed to work with live signals.

%% Beacon Packet Recovery
% The following steps happen sequentially to recover one Non-HT packet:
% 
% * Packet Detection: First a packet must be detected before any processing
% begins. This is accomplished by auto-correlating input symbols. Since the
% front of each 802.11 OFDM packet contains a repetitive structure called
% the L-STF, peaks will occur in the correlation when this packet is
% present. The L-STF field is then extracted and used for coarse frequency
% estimation.
% 
% * Symbol Timing: Once a packet has been detected, future symbols will be
% collected and cross-correlated against to locate the L-LTF. The resulting
% correlation peaks provide accurate timing for OFDM symbol timing. Once
% the full L-LTF is located, it is extracted and used for channel
% estimation, and fine frequency estimation.
%
% * L-SIG Decoding: The first OFDM symbol after the L-LTF is the L-SIG
% field. This field must be recovered and decoded to determine the
% modulation, code rate, and length of the following payload. The
% information is used to first capture the correct amount of data after the
% L-SIG for the complete payload and to decode that information.
% 
% * Payload Decoding: All OFDM symbols after the L-SIG are buffered to a
% length determined by the L-SIG field. After all the symbols have been
% captured they are demodulated and decoded into their source bits. The
% source bits are then evaluated. This evaluation includes CRC checking,
% header and body extraction. If the packet is of subtype beacon, summary
% information such as SSID will be printed for the recovered packet.
% 
% Once a full packet is received or any failures occur during the
% processing chain, the receiver will return to packet detection to search
% for more packets. This process is repeated for the duration of the
% signal.

%% Streaming Process on Captured Data
% In this example an off-the-air capture is processed to recover beacon
% frames. A Wi-Fi signal was captured using an RF interface with one
% receive antenna at a sampling rate of 20 Msps. The captured waveform is
% stored in a compressed baseband file. The file was created using
% <matlab:edit('BasebandWriter.m') BasebandWriter.m>. 
%
% The captured waveform is processed in a streaming fashion. A block of
% samples is pulled in for processing in each iteration. As many valid
% packets are retrieved as possible. <matlab:edit('BasebandReader.m')
% BasebandReader.m> is used to read blocks of samples from the compressed
% baseband file.

% Create an object to stream the data from the file
BR = BasebandReader( ...
    'Filename',       'nonHTBeaconRxData.bb', ...
    'SamplesPerFrame', 80); % Number of samples in 1 OFDM symbol at 20 MHz

%%
% The center frequency, sample rate and number of channels in the captured
% waveform are provided by the BasebandReader object.

disp(['Center frequency: ' num2str(BR.CenterFrequency/1e6) ' MHz']);
disp(['Sample rate: ' num2str(BR.SampleRate/1e6) ' Msps']);
disp(['Number of receive antennas: ' num2str(BR.NumChannels) ...
     sprintf('\n')]);

%%
% A <matlab:edit('nonHTFrontEnd.m') nonHTFrontEnd.m> object performs
% front-end processing and L-SIG decoding. The object is configured with a
% channel bandwidth of 20 MHz to process Non-HT packets. Only one receive
% antenna is supported.
RxFE = nonHTFrontEnd( ....
    'ChannelBandwidth', 'CBW20', ....
    'SymbolTimingThreshold', 0.2);  

%%
% A while loop is used to process blocks of samples and recover beacon
% packets until no more data is available in the baseband file. In each
% iteration of the loop a block of samples are read from the baseband file
% and are processed by |RxFE|. |RxFE| performs front-end processing and
% buffers samples until a packet has been detected and the payload
% received. When |payloadFull| is true, the full payload has been buffered
% and |RxFE| returns variables to allow the data within the packet to be
% recovered:
%
% * |cfgNonHT| is the recovered packet parameters from L-SIG
% * |rxNonHTData| is the time-domain Non-HT data field signal
% * |chanEst| is the channel estimates obtained from the L-LTF
% * |noiseVar| is the fixed noise variance value
%
% The packet payload bits are recovered from the Non-HT data field samples
% using <matlab:doc('wlanNonHTDataRecover') wlanNonHTDataRecover>. The bits
% are then processed by <matlab:edit('nonHTBeaconRxMPDUDecode.m')
% nonHTBeaconRxMPDUDecode.m> to determine whether the recovered packet is a
% valid beacon. If a valid beacon is detected
% <matlab:edit('nonHTBeaconRxOutputDisplay.m')
% nonHTBeaconRxOutputDisplay.m> is used to display its properties.

% A recovery configuration object is used to specify zero-forcing
% equalization for the data recovery
cfgRec = wlanRecoveryConfig('EqualizationMethod', 'ZF');

% Symbol-by-symbol streaming process
numValidPackets = 0;
while ~isDone(BR)
    % Pull in one OFDM symbol, i.e. 80 samples
    data = step(BR);
    
    % Perform front-end processing and payload buffering
    [payloadFull, cfgNonHT, rxNonHTData, chanEst, noiseVar] = ...
        step(RxFE, data);
    
    if payloadFull        
        % Recover payload bits
        recBits = wlanNonHTDataRecover(rxNonHTData, chanEst, ...
            noiseVar, cfgNonHT, cfgRec);
        
        % Evaluate recovered bits
        [validBeacon, MPDU] = nonHTBeaconRxMPDUDecode(recBits);
        if validBeacon
            nonHTBeaconRxOutputDisplay(MPDU); % Display SSID
            numValidPackets = numValidPackets + 1;
        end
    end    
end

disp([num2str(numValidPackets), ' Valid Beacon Packets Found']);

%% Further Exploration
% See <sdruwlanofdm80211BeaconRx.html 802.11(TM) OFDM Beacon Receiver with
% USRP(R) Hardware> for an example of processing live signals with USRP(R).

%% Appendix
% This example uses the following helper functions and objects:
%
% * <matlab:edit('nonHTFrontEnd.m') nonHTFrontEnd.m>
% * <matlab:edit('nonHTBeaconRxMPDUDecode.m') nonHTBeaconRxMPDUDecode.m>
% * <matlab:edit('nonHTBeaconRxOutputDisplay.m') nonHTBeaconRxOutputDisplay.m>
% * <matlab:edit('BasebandReader.m') BasebandReader.m>
% * <matlab:edit('BasebandWriter.m') BasebandWriter.m>

%% Selected Bibliography
% # IEEE Std 802.11(TM)-2012 IEEE Standard for Information technology -
% Telecommunications and information exchange between systems - Local and
% metropolitan area networks - Specific requirements - Part 11: Wireless
% LAN Medium Access Control (MAC) and Physical Layer (PHY) Specifications.
% # IEEE Std 802.11ac(TM)-2013 IEEE Standard for Information technology -
% Telecommunications and information exchange between systems - Local and
% metropolitan area networks - Specific requirements - Part 11: Wireless
% LAN Medium Access Control (MAC) and Physical Layer (PHY) Specifications -
% Amendment 4: Enhancements for Very High Throughput for Operation in Bands
% below 6 GHz.

displayEndOfDemoMessage(mfilename)