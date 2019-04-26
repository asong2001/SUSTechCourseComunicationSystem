function cfgVHTRx = helperVHTConfigRecover(LSIGBits, VHTSIGABits)
% helperVHTConfigRecover Create a configuration object from signaling bits
%
%   CFGVHT = helperVHTConfigRecover(LSIGBITS,VHTSIGABITS) returns a VHT
%   configuration object of type format configuration object of type 
%   <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a> given recovered bits from L-SIG and VHT-SIG-A for a 
%   single user transmission.

% Copyright 2015 The MathWorks, Inc.

cfgVHTRx = wlanVHTConfig;

% Retrieve information from VHT-SIG-A
VHTSIGABits = double(reshape(VHTSIGABits, 24, 2)');
if all(VHTSIGABits(1,1:2) == [0 0])
    cfgVHTRx.ChannelBandwidth = 'CBW20';
    numSD = 52;
elseif  all(VHTSIGABits(1,1:2) == [1 0])
    cfgVHTRx.ChannelBandwidth = 'CBW40';
    numSD = 108;
elseif  all(VHTSIGABits(1,1:2) == [0 1])
    cfgVHTRx.ChannelBandwidth = 'CBW80';
    numSD = 234;
else
    cfgVHTRx.ChannelBandwidth = 'CBW160';
    numSD = 468;
end    
cfgVHTRx.STBC = logical(VHTSIGABits(1, 4)); 
cfgVHTRx.GroupID = bi2de(VHTSIGABits(1, 5:10)); 
cfgVHTRx.NumSpaceTimeStreams = bi2de(VHTSIGABits(1, 11:13)) + 1; 
cfgVHTRx.NumTransmitAntennas = cfgVHTRx.NumSpaceTimeStreams; 
cfgVHTRx.PartialAID = bi2de(VHTSIGABits(1, 14:22)); 
cfgVHTRx.MCS = bi2de(VHTSIGABits(2, 5:8)); 
cfgVHTRx.Beamforming = logical(VHTSIGABits(2, 9));

channelCodingBit = bi2de(VHTSIGABits(2, 3));
if channelCodingBit==1
    % Channel coding is LDPC
    % cfgVHTRx.ChannelCoding = 'LDPC';
    coder.internal.error('wlan:helperVHTConfigRecover:InvalidChannelCoding');
else
    % Channel coding is BCC which is implied by default so do not set
    % cfgVHTRx.ChannelCoding = 'BCC';
end

% Retrieve information from L-SIG
RXTime = (bi2de(double(LSIGBits(6:17)')) + 3)/3*4 + 20; % 4 symbol range, [RXTime - 3: RXTime]

% Derive number of OFDM symbols and APEP/PSDU lengths
numSTS = cfgVHTRx.NumSpaceTimeStreams;
NVHTLTFVec = [1 2 4 4 6 6 8 8]; 
numPreambSym = 9 + NVHTLTFVec(numSTS);    
if VHTSIGABits(2, 1)
    cfgVHTRx.GuardInterval = 'Short';
    numDataSym = floor((RXTime/4 - numPreambSym)*10.0/9.0) - VHTSIGABits(2, 2);
else
    cfgVHTRx.GuardInterval = 'Long';
    numDataSym = RXTime/4- numPreambSym;  % Precise
end
[numDBPS, numES] = getMCSTable(cfgVHTRx.MCS, numSD, ....
    cfgVHTRx.NumSpaceTimeStreams/(1+cfgVHTRx.STBC)); 
numTailBits = 6;
% Calculate received PSDULength and set it to be the APEPLength
cfgVHTRx.APEPLength = floor((numDataSym*numDBPS - numTailBits*numES - 16)/8);

end

function [Ndbps, Nes] = getMCSTable(MCS, Nsd, Nss)

switch MCS
  case 0
    Nbpscs = 1; % 'BPSK'
    rate   = 1/2;
  case 1
    Nbpscs = 2; % 'QPSK'
    rate   = 1/2;
  case 2
    Nbpscs = 2; 
    rate   = 3/4;
  case 3
    Nbpscs = 4; % '16QAM'
    rate   = 1/2;
  case 4
    Nbpscs = 4; 
    rate   = 3/4;
  case 5
    Nbpscs = 6; % '64QAM'
    rate   = 2/3;
  case 6
    Nbpscs = 6; 
    rate   = 3/4;
  case 7
    Nbpscs = 6;
    rate   = 5/6;
  case 8
    Nbpscs = 8; % 256QAM
    rate   = 3/4;
  otherwise % MCS == 9
    Nbpscs = 8;
    rate   = 5/6;
end    

Ndbps = Nsd * Nbpscs * Nss * rate;

% Handle exceptions to Nes generic rule - Table 7.13 [2].
%   For each case listed, work off the Ndbps value and create a look-up
%   table for the Nes value.
%   Only 9360 has a valid value from the generic rule also, 
%   all others are exceptions
NdbpsVec = [2457 8190 9828 9360 14040 9828 16380 19656 21840 14976 22464];
expNes =   [   3    6    6    6     8    6     9    12    12     8    12];

exceptIdx = find(Ndbps == NdbpsVec);
if ~isempty(exceptIdx)
    if (Ndbps == 9360) && (Nss == 5) % One valid case for 160, 80+80
        Nes = 5;
    else  % Two exception cases
        Nes = expNes(exceptIdx(1));
    end
else  % Generic rule: 3.6*600 - for a net 600Mbps per encoder
    Nes = ceil(Ndbps/2160);
end

end

% [EOF]
