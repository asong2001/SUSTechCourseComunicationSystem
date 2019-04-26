%% 802.11ac(TM) Signal Recovery with Preamble Decoding
%
% This example shows how to detect packets and decode payload bits in a
% received IEEE(R) 802.11ac(TM) VHT waveform. The receiver recovers the
% packet format parameters from the preamble fields to decode the data.

% Copyright 2015 The MathWorks, Inc.

%% Introduction
% In a single user 802.11ac packet the transmission parameters are signaled
% to the receiver using the L-SIG and VHT-SIG-A preamble fields [ <#15 1>
% ]:
%
% * The L-SIG field contains information to allow the receiver to determine
% the transmission time of a packet.
% * The VHT-SIG-A field contains the transmission parameters including the
% modulation and coding scheme, number of space-time streams and channel
% coding.
%
% In this example we detect and decode a packet within a generated
% waveform. All transmission parameters apart from the channel bandwidth
% are assumed unknown and are therefore retrieved from the decoded L-SIG
% and VHT-SIG-A preamble fields in each packet. The retrieved transmission
% configuration is used to decode the VHT-SIG-B and VHT Data fields.
% Additionally the following analysis is performed:
%
% * The waveform of the detected packet is displayed.
% * The spectrum of the detected packet is displayed.
% * The constellation of the equalized data symbols per spatial stream are
% displayed.

%% Waveform Transmission
% In this example an 802.11ac VHT single-user waveform is generated locally
% but a captured waveform could be used. MATLAB(R) can be used to acquire
% I/Q data from a wide range of instruments using the Instrument Control
% Toolbox(TM) and software defined radio platforms.
%
% The locally generated waveform is impaired by a 3x3 TGac fading channel,
% additive white Gaussian noise and carrier frequency offset. To generate a
% waveform locally we configure a VHT packet format configuration object.
% Note that the VHT packet configuration object is used at the transmitter
% side only. The receiver will dynamically formulate another VHT
% configuration object when the packet is decoded. The helper function
% <matlab:edit('vhtSigRecGenerateWaveform.m') vhtSigRecGenerateWaveform>
% generates the impaired waveform locally. The processing steps within the
% helper function are:
% 
% * A PSDU is randomly generated and encoded into a VHT waveform
% * The waveform is passed through a TGac fading channel model
% * Carrier frequency offset is added to the waveform
% * Additive white Gaussian noise is added to the waveform

% VHT link parameters
cfgVHTTx = wlanVHTConfig( ...
    'ChannelBandwidth',    'CBW80', ...
    'NumTransmitAntennas', 3, ...
    'NumSpaceTimeStreams', 2, ...
    'SpatialMapping',      'Hadamard', ...
    'STBC',                true, ...
    'MCS',                 5, ...
    'GuardInterval',       'Long', ...
    'APEPLength',          1050);

% Propagation channel
numRx = 3;                  % Number of receive antennas
delayProfile = 'Model-C';   % TGac channel delay profile

% Impairments
noisePower = -30;  % Noise power to apply in dBW
cfo = 62e3;        % Carrier frequency offset (Hz)

% Generated waveform parameters
numTxPkt = 1;      % Number of transmitted packets
idleTime = 20e-6;  % Idle time before and after each packet

% Generate waveform
rxWave = vhtSigRecGenerateWaveform(cfgVHTTx, numRx, ...
    delayProfile, noisePower, cfo, numTxPkt, idleTime);

%% Packet Recovery
% The signal to process is stored in the variable |rxWave|. The processing
% steps to recover a packet are:
% 
% * The packet is detected and synchronized
% * The L-SIG field is extracted and its information bits are recovered to
% determine the length of the packet in microseconds
% * The VHT-SIG-A field is extracted and its information bits are recovered
% * The packet format parameters are retrieved from the decoded L-SIG and
% VHT-SIG-A bits
% * The VHT-LTF field is extracted to perform MIMO channel estimation for
% decoding the VHT-SIG-B and VHT Data fields
% * The VHT-SIG-B field is extracted and its information bits recovered
% * The VHT-Data field is extracted and the PSDU and VHT-SIG-B CRC bits
% recovered using the retrieved packet parameters
%
% The start and end indices for some preamble fields are independent of
% transmission parameters excluding the channel bandwidth. These indices
% are calculated using a default transmission configuration object with the
% known bandwidth.

cfgVHTRx = wlanVHTConfig('ChannelBandwidth', cfgVHTTx.ChannelBandwidth);
idxLSTF = wlanFieldIndices(cfgVHTRx, 'L-STF'); 
idxLLTF = wlanFieldIndices(cfgVHTRx, 'L-LTF'); 
idxLSIG = wlanFieldIndices(cfgVHTRx, 'L-SIG'); 
idxSIGA = wlanFieldIndices(cfgVHTRx, 'VHT-SIG-A'); 

%%
% The following code configures objects and variables for processing.

chanBW = cfgVHTTx.ChannelBandwidth;
sr = helperSampleRate(chanBW);

% Setup plots for example
[SpectrumAnalyzer, TimeScope, ConstellationDiagram] = ...
    vhtSigRecSetupPlots(sr);

% Minimum packet length is 10 OFDM symbols
lstfLen = double(idxLSTF(2)); % Number of samples in L-STF
minPktLen = lstfLen*5;

rxWaveLen = size(rxWave, 1);

%% Front-End Processing
% The front-end processing consists of packet detection, coarse carrier
% frequency offset correction, timing synchronization and fine carrier
% frequency offset correction. A <matlab:doc('while') while> loop is used
% to detect and synchronize a packet within the received waveform. The
% sample offset |searchOffset| is used to index into |rxWave| to detect a
% packet. The first packet within |rxWave| is detected and processed. If
% the synchronization fails for the detected packet, the sample index
% offset |searchOffset| is incremented to move beyond the processed packet
% in |rxWave|. This is repeated until a packet has been successfully
% detected and synchronized.

searchOffset = 0; % Offset from start of waveform in samples
while (searchOffset + minPktLen) <= rxWaveLen
    % Packet detection
    pktStartOffset = helperPacketDetect(rxWave(1+searchOffset:end, :), ...
        chanBW) - 1;

    % Adjust packet offset
    pktStartOffset = searchOffset + pktStartOffset;
    if isempty(pktStartOffset) || (pktStartOffset + idxLSIG(2) > rxWaveLen)
        error('** No packet detected **');
    end

    % Coarse frequency offset estimation using L-STF
    LSTF = rxWave(pktStartOffset + (idxLSTF(1):idxLSTF(2)), :);
    coarseFreqOffset = wlanCoarseCFOEstimate(LSTF, chanBW);

    % Coarse frequency offset compensation
    rxWave(pktStartOffset+1:end,:) = helperFrequencyOffset( ...
        rxWave(pktStartOffset+1:end,:), sr, -coarseFreqOffset);

    % Symbol timing synchronization: 4 OFDM symbols to search for L-LTF
    LLTFSearchBuffer = rxWave(pktStartOffset + idxLSTF(2)/2 + ...
        (idxLSTF(1):idxLLTF(2)), :);    
    LLTFStartOffset = helperSymbolTiming(LLTFSearchBuffer, chanBW) - 1;

    % If no L-LTF detected skip samples and continue searching
    if isempty(LLTFStartOffset)
        fprintf('** No L-LTF detected **\n\n');
        searchOffset = pktStartOffset+lstfLen; 
        continue; 
    end

    % Adjust the packet offset now that the L-LTF offset is known
    pktStartOffset = pktStartOffset + LLTFStartOffset - ...
        double(idxLSTF(2)/2);
    if (pktStartOffset < 0) || ((pktStartOffset + minPktLen) > rxWaveLen) 
        fprintf('** Timing offset invalid **\n\n');
        searchOffset = pktStartOffset+lstfLen; 
        continue; 
    end
    
    % Timing synchronization complete: packet detected
    fprintf('Packet detected at index %d\n\n',pktStartOffset+1);
    
    % Fine frequency offset estimation using L-LTF
    LLTF = rxWave(pktStartOffset + (idxLLTF(1):idxLLTF(2)), :);
    fineFreqOffset = wlanFineCFOEstimate(LLTF, chanBW);

    % Fine frequency offset compensation
    rxWave(pktStartOffset+1:end,:) = helperFrequencyOffset( ...
    rxWave(pktStartOffset+1:end,:), sr, -fineFreqOffset);

    % Display estimated carrier frequency offset
    cfoCorrection = coarseFreqOffset + fineFreqOffset; % Total CFO
    fprintf('Estimated CFO: %5.1f Hz\n\n',cfoCorrection);
    
    break; % Front-end processing complete, stop searching for a packet
end

%% L-SIG Decoding
% In a VHT transmission the L-SIG field is used to determine the receive
% time, or RXTIME, of the packet. RXTIME is calculated using the field bits
% of the L-SIG payload [ <#15 1> Eq. 22-105]. The number of samples which
% contain the packet within |rxWave| can then be calculated. The L-SIG
% payload is decoded using an estimate of the channel and noise power
% obtained from the L-LTF.

% Channel estimation using L-LTF
LLTF = rxWave(pktStartOffset + (idxLLTF(1):idxLLTF(2)), :);
demodLLTF = wlanLLTFDemodulate(LLTF, chanBW);
chanEstLLTF = wlanLLTFChannelEstimate(demodLLTF, chanBW);

% Estimate noise power in NonHT fields
noiseVarNonHT = helperNoiseEstimate(demodLLTF);

% Recover L-SIG field bits
disp('Decoding L-SIG... ');
[recLSIGBits, failCheck] = wlanLSIGRecover( ...
    rxWave(pktStartOffset + (idxLSIG(1):idxLSIG(2)), :), ...
    chanEstLLTF, noiseVarNonHT, chanBW);

if failCheck % Skip L-STF length of samples and continue searching
    disp('** L-SIG check fail **');
else
    disp('L-SIG check pass');
end

% Calculate the receive time and corresponding number of samples in the
% packet
lengthBits = recLSIGBits(6:17).';
RXTime = ceil((bi2de(double(lengthBits)) + 3)/3) * 4 + 20; % us
numRxSamples = RXTime * 1e-6 * sr; % Number of samples in receive time

fprintf('RXTIME: %dus\n', RXTime);
fprintf('Number of samples in packet: %d\n\n', numRxSamples);

%%
% The waveform and spectrum of the detected packet within |rxWave| are
% displayed given the calculated RXTIME and corresponding number of
% samples.

sampleOffset = max((-lstfLen + pktStartOffset), 1); % First index to plot
sampleSpan = numRxSamples + 2*lstfLen;           % Number samples to plot
% Plot as much of the packet (and extra samples) as we can
plotIdx = sampleOffset:min(sampleOffset + sampleSpan, rxWaveLen);

% Configure TimeScope to display the packet
TimeScope.TimeSpan = sampleSpan/sr;
TimeScope.TimeDisplayOffset = sampleOffset/sr;
TimeScope.YLimits = [0 max(abs(rxWave(:)))];
step(TimeScope, abs(rxWave(plotIdx ,:)));

% Display the spectrum of the detected packet
step(SpectrumAnalyzer, rxWave(pktStartOffset + (1:numRxSamples), :));

%% VHT-SIG-A Decoding
% The VHT-SIG-A field contains the transmission configuration of the
% packet. The VHT-SIG-A bits are recovered using the channel and noise
% power estimates obtained from the L-LTF.

% Recover VHT-SIG-A field bits
disp('Decoding VHT-SIG-A... ');
[recVHTSIGABits, failCRC] = wlanVHTSIGARecover( ...
    rxWave(pktStartOffset + (idxSIGA(1):idxSIGA(2)), :), ...
    chanEstLLTF, noiseVarNonHT, chanBW);

if failCRC
    disp('** VHT-SIG-A CRC fail **');
else
    disp('VHT-SIG-A CRC pass');
end

%% 
% The helper function <matlab:edit('helperVHTConfigRecover.m')
% helperVHTConfigRecover.m> returns a VHT format configuration object,
% |cfgVHTRx|, based on recovered VHT-SIG-A and L-SIG bits. Properties which
% are not required to decode the waveform are set to default values for a
% <matlab:doc('wlanVHTConfig') wlanVHTConfig> object and therefore may
% differ from the value in |cfgVHTTx|. Examples of such properties include
% |NumTransmitAntennas| and |SpatialMapping|.

% Create a VHT format configuration object by retrieving packet parameters
% from the decoded L-SIG and VHT-SIG-A bits
cfgVHTRx = helperVHTConfigRecover(recLSIGBits, recVHTSIGABits);

% Display the transmission configuration obtained from VHT-SIG-A
vhtSigRecDisplaySIGAInfo(cfgVHTRx);

%%
% The information provided by VHT-SIG-A allows the location of subsequent
% fields within the received waveform to be calculated.

% Obtain starting and ending indices for VHT-LTF and VHT-Data fields
% using retrieved packet parameters
idxVHTLTF  = wlanFieldIndices(cfgVHTRx, 'VHT-LTF');
idxVHTSIGB = wlanFieldIndices(cfgVHTRx, 'VHT-SIG-B');
idxVHTData = wlanFieldIndices(cfgVHTRx, 'VHT-Data');

% Warn if waveform does not contain whole packet
if (pktStartOffset + double(idxVHTData(2))) > rxWaveLen 
    fprintf('** Not enough samples to recover entire packet **\n\n');
end

%% VHT-SIG-B Decoding
% The primary use of VHT-SIG-B is for signaling user information in a
% multi-user packet. In a single user packet the VHT-SIG-B carries the
% length of the packet which can also be calculated using the L-SIG and
% VHT-SIG-A (which is demonstrated in the sections above). Despite not
% being required to decode a single user packet, the VHT-SIG-B is recovered
% below and the bits interpreted. The VHT-SIG-B symbols are demodulated
% using a MIMO channel estimate obtained from the VHT-LTF. Note the CRC for
% VHT-SIG-B is carried in the VHT Data field.

% Estimate MIMO channel using VHT-LTF and retrieved packet parameters
demodVHTLTF = wlanVHTLTFDemodulate( ...
    rxWave(pktStartOffset + (idxVHTLTF(1):idxVHTLTF(2)), :), cfgVHTRx);
chanEstVHTLTF = wlanVHTLTFChannelEstimate(demodVHTLTF, cfgVHTRx);

% Estimate noise power in VHT fields
noiseVarVHT = helperNoiseEstimate(demodLLTF, cfgVHTRx);

% VHT-SIG-B Recover
disp('Decoding VHT-SIG-B...');
[sigbBits, sigbSym] = wlanVHTSIGBRecover( ...
    rxWave(pktStartOffset + (idxVHTSIGB(1):idxVHTSIGB(2)),:), ...
        chanEstVHTLTF, noiseVarVHT, chanBW);

% Interpret VHT-SIG-B bits to recover the APEP length (rounded up to a
% multiple of four bytes) and generate reference CRC bits
[sigbAPEPLength, refSIGBCRC] = ...
    vhtSigRecInterpretVHTSIGBBits(sigbBits, cfgVHTRx);
    disp('Decoded VHT-SIG-B contents: ');
    fprintf('  APEP Length (rounded up to 4 byte multiple): %d bytes\n\n', ...
        sigbAPEPLength);

%% VHT Data Decoding
% The reconstructed VHT configuration object can then be used to recover
% the VHT Data field. This includes the VHT-SIG-B CRC bits and PSDU.
%
% The recovered VHT data symbols can then be analyzed as required. In this
% example the equalized constellation of the recovered VHT data symbols per
% spatial stream are displayed.

% Recover PSDU bits using retrieved packet parameters and channel
% estimates from VHT-LTF
disp('Decoding VHT Data field...');
[rxPSDU, rxSIGBCRC, eqDataSym] = wlanVHTDataRecover( ...
    rxWave(pktStartOffset + (idxVHTData(1):idxVHTData(2)), :), ...
    chanEstVHTLTF, noiseVarVHT, cfgVHTRx);

% Plot equalized constellation for each spatial stream
refConst = helperConstellationSymbols(cfgVHTRx);
[Nsd, Nsym, Nss] = size(eqDataSym);
eqDataSymPerSS = reshape(eqDataSym, Nsd*Nsym, Nss);
for iss = 1:Nss
    ConstellationDiagram{iss}.ReferenceConstellation = refConst;
    step(ConstellationDiagram{iss}, eqDataSymPerSS(:, iss));
end

%%
% The CRC bits for VHT-SIG-B recovered in VHT Data are then compared to the
% locally generated reference to determine whether the VHT-SIG-B and VHT
% data service bits have been recovered successfully.

% Test VHT-SIG-B CRC from service bits within VHT Data against
% reference calculated with VHT-SIG-B bits
if ~isequal(refSIGBCRC, rxSIGBCRC)
    disp('** VHT-SIG-B CRC fail **');
else
    disp('VHT-SIG-B CRC pass');
end

%% Appendix
% This example uses the following helper functions:
%
% * <matlab:edit('helperConstellationSymbols.m') helperConstellationSymbols.m>
% * <matlab:edit('helperNoiseEstimate.m') helperNoiseEstimate.m>
% * <matlab:edit('helperPacketDetect.m') helperPacketDetect.m>
% * <matlab:edit('helperSampleRate.m') helperSampleRate.m>
% * <matlab:edit('helperSymbolTiming.m') helperSymbolTiming.m>
% * <matlab:edit('helperVHTConfigRecover.m') helperVHTConfigRecover.m>
% * <matlab:edit('vhtSigRecDisplaySIGAInfo.m') vhtSigRecDisplaySIGAInfo.m>
% * <matlab:edit('vhtSigRecGenerateWaveform.m') vhtSigRecGenerateWaveform.m>
% * <matlab:edit('vhtSigRecInterpretVHTSIGBBits.m') vhtSigRecInterpretVHTSIGBBits.m>
% * <matlab:edit('vhtSigRecSetupPlots.m') vhtSigRecSetupPlots.m>

%% Selected Bibliography
% # IEEE Std 802.11ac(TM)-2013 IEEE Standard for Information technology -
% Telecommunications and information exchange between systems - Local and
% metropolitan area networks - Specific requirements - Part 11: Wireless
% LAN Medium Access Control (MAC) and Physical Layer (PHY) Specifications -
% Amendment 4: Enhancements for Very High Throughput for Operation in Bands
% below 6 GHz.

displayEndOfDemoMessage(mfilename)