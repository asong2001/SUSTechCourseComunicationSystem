function [y, bits] = wlanVHTSIGA(cfgVHT)
%WLANVHTSIGA VHT Signal A (VHT-SIG-A) field
%
%   [Y, BITS] = wlanVHTSIGA(CFGVHT) generates the VHT Signal A
%   (VHT-SIG-A) field time-domain waveform for the VHT transmission format.
%
%   Y is the time-domain VHT-SIG-A field signal. It is a complex matrix of
%   size Ns-by-Nt, where Ns represents the number of time-domain samples
%   and Nt represents the number of transmit antennas.
%
%   BITS is the VHT-SIG-A signaling bits. It is an int8-typed, binary column
%   vector of length 48.
%
%   CFGVHT is the format configuration object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a> which
%   specifies the parameters for the VHT format.
%
%   Example:
%   % Generate the VHT-SIG-A waveform for a VHT 80 MHz transmission format
%
%     cfgVHT = wlanVHTConfig;               % VHT Format configuration
%     cfgVHT.ChannelBandwidth = 'CBW80';    % Set to 80 MHz
%     vSigAOut = wlanVHTSIGA(cfgVHT);
%     %  returns a complex output of 640 samples for the two OFDM symbols
%     %  of length 320 samples each.
%
%   See also wlanVHTConfig, wlanLSIG, wlanVHTSTF, wlanVHTSIGARecover.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

validateattributes(cfgVHT, {'wlanVHTConfig'}, {'scalar'}, mfilename, ...
                   'VHT format configuration object');
cfgInfo = validateConfig(cfgVHT, 'SpatialMCSGID');

% VHT-SIG-A1 structure - Table 22-12, IEEE Std 802.11ac-2013
% BW field
switch cfgVHT.ChannelBandwidth
    case 'CBW20'
        bw = [0 0];
    case 'CBW40'
        bw = [1 0];     % right-msb orientation
    case 'CBW80'
        bw = [0 1];     % right-msb orientation
    otherwise % 'CBW160', 'CBW80+80'
        bw = [1 1];
end

% STBC, NSTS/Partial AID fields
if cfgVHT.NumUsers == 1
    STBC = double(cfgVHT.STBC);
    STSAndPAID = [de2bi(cfgVHT.NumSpaceTimeStreams(1)-1, 3), ...
                  de2bi(cfgVHT.PartialAID, 9)];
else
    STBC = 0;
    STSAndPAID = zeros(1, 12);
    for u = 1:cfgVHT.NumUsers
        STSAndPAID(3*cfgVHT.UserPositions(u)+(1:3)) = ...
                    de2bi(cfgVHT.NumSpaceTimeStreams(u), 3);
    end
end

% Preset TransmitPowerSaveNotAllowed to false
TransmitPowerSaveNotAllowed = 0;

% Assemble fields with reserved bits
vhtsiga1 = int8([bw 1 STBC de2bi(cfgVHT.GroupID, 6), ...
            STSAndPAID, TransmitPowerSaveNotAllowed 1].');

% VHT-SIG-A2 structure - Table 22-12, IEEE Std 802.11ac-2013
% Guard interval bits
if strcmp(cfgVHT.GuardInterval, 'Long')
    b0_2  = 0;
    b1_2  = 0;
else % Short GI
    Nsym = cfgInfo.NumDataSymbols;
    b0_2 = 1;
    b1_2 = double((mod(Nsym, 10)==9));     
end

% Channel coding bits
% Update b3 for LDPC coding, once enabled
if cfgVHT.NumUsers == 1
    if strcmp(cfgVHT.ChannelCoding, 'BCC')
        b2_2 = 0;
        b3_2 = 0;
    else  % LDPC
        b2_2 = 1;
        b3_2 = 0;     
    end
    b2Tob7 = [b2_2; b3_2; de2bi(cfgVHT.MCS(1), 4, 'right-msb').'];
else 
    MUCoding = ones(4, 1);
    for u = 1:cfgVHT.NumUsers
        MUCoding(cfgVHT.UserPositions(u)+1) = ...
                        double(strcmp(cfgVHT.ChannelCoding, 'LDPC'));
    end
    b2Tob7 = [MUCoding(1); 0; MUCoding(2:4); 1];
end

% Set BEAMFORMED bit
b8 = ((cfgVHT.NumUsers == 1) && strcmp(cfgVHT.SpatialMapping, 'Custom') ...
      && cfgVHT.Beamforming) || (cfgVHT.NumUsers > 1);

% Concatenate the first 0-9 bits
vhtsiga2_09 = int8([b0_2; b1_2; b2Tob7; b8; 1]);

% Generate the CRC
crc = wlanCRCGenerate([vhtsiga1; vhtsiga2_09]);

% VHT-SIG-A2 bits
vhtsiga2 = [vhtsiga2_09; crc; zeros(6,1,'int8')]; % 24 bits

% Concatenate the SIG-A1 and A2 fields together - 48 bits
bits = [vhtsiga1; vhtsiga2];

%% Process VHT-SIG-A bits

% Get OFDM parameters
cfgOFDM = wlanGetOFDMConfig(cfgVHT.ChannelBandwidth, 'Long', 'Legacy');
FFTLen  = cfgOFDM.FFTLength;
CPLen   = cfgOFDM.CyclicPrefixLength;
num20   = FFTLen/64;

encodedSIGA  = wlanBCCEncode(bits, '1/2');
interleavedSIGA1 = wlanBCCInterleave(encodedSIGA(1:48), 'NON_HT', 48, 1);
interleavedSIGA2 = wlanBCCInterleave(encodedSIGA(49:end), 'NON_HT', 48, 1);
firstDataSym = wlanConstellationMapper(interleavedSIGA1, 1);
secDataSym   = wlanConstellationMapper(interleavedSIGA2, 1, pi/2);

firstSym = complex(zeros(FFTLen, 1)); 
secSym   = complex(zeros(FFTLen, 1)); 
% Add pilot subcarriers, IEEE Std 802.11ac-2013, Eqn 22-28
firstSym(cfgOFDM.DataIndices)  = repmat(firstDataSym, num20, 1);
Nsym = 1; % Create pilots symbol-by-symbol
z = 1; % Offset by 1 to account for L-SIG pilot symbol
firstSym(cfgOFDM.PilotIndices) = repmat(nonHTPilots(Nsym, z), num20, 1);
secSym(cfgOFDM.DataIndices)    = repmat(secDataSym, num20, 1);
z = 2; % Offset by 2 to account for L-SIG and first VHT-SIG-A pilot symbols
secSym(cfgOFDM.PilotIndices)   = repmat(nonHTPilots(Nsym, z), num20, 1);

% Replicate over bandwidth, with tone rotation (gamma)
firstSymAll = firstSym .* cfgOFDM.CarrierRotations;
secSymAll   = secSym   .* cfgOFDM.CarrierRotations;

% Concatenate the two symbols
SymAll = [firstSymAll, secSymAll];

% Cyclic shift addition
% The format is set to OFDM due to legacy mode. The cyclic shift is
% applied on each transmit antenna.
numTx = cfgVHT.NumTransmitAntennas;
vhtCycShift = complex(zeros(FFTLen, 2, numTx));
csh = getCyclicShiftVal('OFDM', numTx, 20*num20);
for i = 1:2
    % Replicate VHT-SIG-A field over multiple antennas
    vhtsigMIMO = repmat(SymAll(:,i), 1, numTx);
    
    vhtCycShift(:,i,:) = wlanCyclicShift(vhtsigMIMO, csh, FFTLen, 'Tx');
end

% OFDM modulate
wout = wlanOFDMModulate(vhtCycShift, CPLen);
y  = wout * cfgOFDM.NormalizationFactor / sqrt(numTx);

end