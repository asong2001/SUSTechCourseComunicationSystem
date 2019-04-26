function [y, bits] = wlanHTSIG(cfgHT)
%WLANHTSIG HT SIGNAL (HT-SIG) field
%
%   [Y, BITS] = wlanHTSIG(CFGHT) generates the HT SIGNAL (HT-SIG)
%   field time-domain waveform for the HT-Mixed transmission format.
%
%   Y is the time-domain HT-SIG field signal. It is a complex matrix of
%   size Ns-by-Nt, where Ns represents the number of time-domain samples
%   and Nt represents the number of transmit antennas.
%
%   BITS is the signaling bits used for the HT SIGNAL field. It is an
%   int8-typed, binary column vector of length 48.
%
%   CFGHT is the format configuration object of type <a href="matlab:help('wlanHTConfig')">wlanHTConfig</a> which
%   specifies the parameters for the HT-Mixed format.
%
%   Example: 
%   %  Generate the HT-SIG waveform for a HT 40MHz transmission format
%
%      cfgHT = wlanHTConfig;                % Format configuration
%      cfgHT.ChannelBandwidth = 'CBW40';    % Set to 40MHz
%      htSigOut = wlanHTSIG(cfgHT);
%
%   See also wlanHTConfig, wlanLSIG, wlanHTSTF, wlanHTSIGRecover.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

% Validate input HT format
validateattributes(cfgHT, {'wlanHTConfig'}, {'scalar'}, mfilename, ...
                   'HT-Mixed format configuration object');
validateConfig(cfgHT, 'MCSSTSTx'); % MCS, Length, STS, Tx              

%% Build the signaling bits
% HT-SIG1 structure
b17 = de2bi(cfgHT.MCS, 7, 'right-msb').';
switch cfgHT.ChannelBandwidth
    case {'CBW40'}
        b8 = 1;
    otherwise
        b8 = 0; % for 20MHz only, offset not implemented
end
b924 = de2bi(cfgHT.PSDULength, 16, 'right-msb').';

htsig1 = [b17; b8; b924]; % 24 bits

% HT-SIG2 structure
b1 = double(cfgHT.RecommendSmoothing);
b2 = double(cfgHT.PSDULength~=0);    % Not Sounding

% Set Aggregation bit to always be false (0)
b4 = 0;

% STBC value
Nss = floor(cfgHT.MCS/8)+1;
STBC = cfgHT.NumSpaceTimeStreams - Nss;
b56 = de2bi(STBC, 2, 'right-msb').';

b7 = double(strcmp(cfgHT.ChannelCoding, 'LDPC'));
b8 = double(strcmp(cfgHT.GuardInterval, 'Short'));
if inESSMode(cfgHT)
    numESS = cfgHT.NumExtensionStreams;
else
    numESS = 0;
end
b910 = de2bi(numESS, 2, 'right-msb').';

% Concatenate the first 0-9 bits
htsig2_09 = [b1; b2; 1; b4; b56; b7; b8; b910];

% Generate the CRC
crc = wlanCRCGenerate([htsig1; htsig2_09]);

% HT SIG2 bits
htsig2 = [htsig2_09; crc; zeros(6,1, 'int8')]; % 24 bits

% Concatenate the HT-SIG-1 and HT-SIG-2 fields together: 48 bits
bits = [htsig1; htsig2];

%% Process HT-SIG bits
encodedSIG  = wlanBCCEncode(bits, '1/2');
%   Update Interleave to accept matrices as input
interleavedSIG1 = wlanBCCInterleave(encodedSIG(1:48), 'NON_HT', 48, 1);
interleavedSIG2 = wlanBCCInterleave(encodedSIG(49:end), 'NON_HT', 48, 1);
dataSym = wlanConstellationMapper([interleavedSIG1; interleavedSIG2], ...
                                  1, pi/2);

% Get OFDM parameters
cfgOFDM = wlanGetOFDMConfig(cfgHT.ChannelBandwidth, 'Long', 'Legacy');
FFTLen  = cfgOFDM.FFTLength;
CPLen   = cfgOFDM.CyclicPrefixLength;
num20   = FFTLen/64;

firstSym = complex(zeros(FFTLen,1)); 
secSym   = complex(zeros(FFTLen,1)); 

% Add pilot subcarriers, from IEEE Std 802.11-2012, Eqn 20-17
firstSym(cfgOFDM.DataIndices)   = repmat(dataSym(1:48), num20, 1);
Nsym = 1; % Create pilots symbol-by-symbol
z = 1;    % Offset by 1 to account for L-SIG pilot symbol
firstSym(cfgOFDM.PilotIndices) = repmat(nonHTPilots(Nsym, z), num20, 1);
secSym(cfgOFDM.DataIndices)    = repmat(dataSym(49:end), num20, 1);
z = 2; % Offset by 2 to account for L-SIG and first HT-SIG pilot symbols
secSym(cfgOFDM.PilotIndices)   = repmat(nonHTPilots(Nsym, z), num20, 1);

% Replicate over bandwidth, with tone rotation (gamma)
firstSymAll = firstSym .* cfgOFDM.CarrierRotations;
secSymAll   = secSym   .* cfgOFDM.CarrierRotations;

% Concatenate the two HT-SIG symbols
SymAll = [firstSymAll, secSymAll];

% Cyclic shift addition
numTx = cfgHT.NumTransmitAntennas;
htCycShift = complex(zeros(FFTLen, 2, numTx));
%   The FORMAT is set to OFDM due to legacy mode. 
csh = getCyclicShiftVal('OFDM', numTx, 20*num20);
for i = 1:2
    % Replicate HTSIG field over multiple transmit antennas
    htsigMIMO = repmat(SymAll(:,i), 1, cfgHT.NumTransmitAntennas);    
    htCycShift(:,i,:) = wlanCyclicShift(htsigMIMO, csh, FFTLen, 'Tx');
end

wout = wlanOFDMModulate(htCycShift, CPLen);
y  = wout * cfgOFDM.NormalizationFactor / sqrt(numTx);

% [EOF]
