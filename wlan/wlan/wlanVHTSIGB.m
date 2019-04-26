function [y, bits]= wlanVHTSIGB(cfgVHT)
%WLANVHTSIGB VHT Signal B (VHT-SIG-B) field
%
%   [Y, BITS] = wlanVHTSIGB(CFGVHT) generates the VHT Signal B (VHT-SIG-B)
%   field time-domain waveform for the VHT transmission format.
%
%   Y is the time-domain VHT-SIG-B field signal. It is a complex matrix of
%   size Ns-by-Nt where Ns represents the number of time-domain samples
%   and Nt represents the number of transmit antennas.
%
%   BITS is the non-repeated signaling bits used for the VHT-SIG-B field.
%   It is an int8-typed, real matrix of size Nb-by-Nu, where Nb is 26 for
%   20 MHz, 27 for 40 MHz, and 29 for 80 MHz and 160 MHz channel
%   bandwidths, respectively, and Nu is the number of users.
%
%   CFGVHT is the format configuration object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a> which
%   specifies the parameters for the VHT format.
%
%   Example:
%   % Generate the VHT-SIG-B waveform for a VHT 80 MHz transmission format
%
%     cfgVHT = wlanVHTConfig;                % Format configuration
%     cfgVHT.ChannelBandwidth = 'CBW80';     % Set to 80 MHz
%     vSigBOut = wlanVHTSIGB(cfgVHT);
%
%   See also wlanVHTConfig, wlanVHTLTF, wlanVHTData, wlanVHTSIGBRecover.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

validateattributes(cfgVHT, {'wlanVHTConfig'}, {'scalar'}, mfilename, ...
                   'VHT format configuration object');
validateConfig(cfgVHT, 'SMappingMCS');

chanBW      = cfgVHT.ChannelBandwidth;
numUsers    = cfgVHT.NumUsers; 
numSTS      = cfgVHT.NumSpaceTimeStreams;
numSTSTotal = sum(numSTS);

% Set up constants related to channel bandwidth
%   Table 22-5, 22-14, IEEE Std 802.11ac-2013
switch chanBW
  case 'CBW20'
    numSD = 52;
    numSIGBBits = 26;
  case 'CBW40'
    numSD = 108;
    numSIGBBits = 27;
  case 'CBW80'
    numSD = 234;
    numSIGBBits = 29;
  otherwise  % 'CBW160'
    numSD = 468;
    numSIGBBits = 29;
end

% Fill data subcarriers for all users
APEPLen = repmat(cfgVHT.APEPLength, 1, numUsers/length(cfgVHT.APEPLength));            
vecMCS  = repmat(cfgVHT.MCS, 1, numUsers/length(cfgVHT.MCS));
data    = complex(zeros(numSD, numSTSTotal));
bits    = zeros(numSIGBBits, numUsers, 'int8');

for u = 1:cfgVHT.NumUsers
    [dataForThisUser, bits(:, u)] = getDataSubcarrierPerUser(chanBW, ...
                                numUsers, APEPLen(u), vecMCS(u));
    data(:, sum(numSTS(1:u-1))+(1:numSTS(u))) = repmat(dataForThisUser, ...
                                                       1, numSTS(u));
end

% Apply the first column of the PVHTLTF matrix
if any(numSTSTotal == [4 7 8])
    % Flip the 4th and 8th STS
    % For all other numSTS, PVHTLTF first column is all ones.
    data(:,4:4:end) = -data(:,4:4:end); 
end

% Get OFDM parameters
cfgOFDM = wlanGetOFDMConfig(chanBW, 'Long', 'VHT', numSTSTotal);
FFTLen = cfgOFDM.FFTLength;
CPLen  = cfgOFDM.CyclicPrefixLength;

% Tone packing for different CBW for a single stream, single symbol
symbol = complex(zeros(FFTLen, numSTSTotal));
symbol(cfgOFDM.DataIndices, :) = data;

% Generate pilot sequence, from Eqn 22-47, IEEE Std 802.11ac-2013
Nsym = 1; % One OFDM symbol in VHT-SIG-B
z = 3; % Offset by 3 to allow for L-SIG and VHT-SIG-A pilot symbols
pilots = vhtPilots(Nsym,z,chanBW,numSTSTotal); % Pilots: Nsp-by-1-by-Nsts
symbol(cfgOFDM.PilotIndices, :) = pilots(:,1,:);

% Tone rotation (gamma) - based on legacy Nfft
symbol = bsxfun(@times, symbol, cfgOFDM.CarrierRotations);

% Cyclic shift addition
csh = getCyclicShiftVal('VHT', numSTSTotal, FFTLen/64*20);
vhtsigCycShift = wlanCyclicShift(symbol, csh, FFTLen, 'Tx');

vhtsigSpatialMapped = wlanSpatialMapping(vhtsigCycShift, ...
    cfgVHT.SpatialMapping, cfgVHT.NumTransmitAntennas, ...
    cfgVHT.SpatialMappingMatrix);

wout = wlanOFDMModulate(reshape(vhtsigSpatialMapped, FFTLen, 1, ...
                                cfgVHT.NumTransmitAntennas), CPLen);

y = wout * cfgOFDM.NormalizationFactor;

end

%--------------------------------------------------------------------------
function [dataPerUser, bitsPerUser] = getDataSubcarrierPerUser(chanBW, ...
                                            numUsers, APEPLength, MCS)
%

% Set up the bits
if APEPLength == 0 % NDP support for a single user
    % IEEE Std 802.11ac-2013, Table 22-15.
    switch chanBW
        case 'CBW20' % 20
            sigbBits = [0 0 0 0 0 1 1 1 0 1 0 0 0 1 0 0 0 0 1 0];
        case 'CBW40' % 21
            sigbBits = [1 0 1 0 0 1 0 1 1 0 1 0 0 0 1 0 0 0 0 1 1];
        otherwise % {'CBW80', 'CBW80+80', 'CBW160'} % 23
            sigbBits = [0 1 0 1 0 0 1 1 0 0 1 0 1 1 1 1 1 1 1 0 0 1 0];
    end
    
    bitsPerUser = [sigbBits zeros(1,6,'int8')].';
else
    
    % Bit values as per IEEE Std 802.11ac-2013, Table 22-14. 
    %   Right-msb orientation for the length bits
    APEPLenOver4 = ceil(APEPLength/4);
    if 1 == numUsers % SU PPDU allocation
        switch chanBW
            case 'CBW20' % 26
                bitsPerUser = [de2bi(APEPLenOver4, 17) ones(1,3, 'int8') ...
                    zeros(1,6,'int8')].';
            case 'CBW40' % 27
                bitsPerUser = [de2bi(APEPLenOver4, 19) ones(1,2, 'int8') ...
                    zeros(1,6,'int8')].';
            otherwise    % 29 for {'CBW80', 'CBW80+80', 'CBW160'}
                bitsPerUser = [de2bi(APEPLenOver4, 21) ones(1,2, 'int8') ...
                    zeros(1,6,'int8')].';
        end
    else % MU PPDU allocation
        switch chanBW
            case 'CBW20' % 26
                bitsPerUser = [de2bi(APEPLenOver4, 16) de2bi(MCS, 4) ...
                    zeros(1,6,'int8')].';
            case 'CBW40' % 27
                bitsPerUser = [de2bi(APEPLenOver4, 17) de2bi(MCS, 4) ...
                    zeros(1,6,'int8')].';
            otherwise    % 29 for {'CBW80', 'CBW80+80', 'CBW160'}
                bitsPerUser = [de2bi(APEPLenOver4, 19), de2bi(MCS, 4) ...
                    zeros(1,6,'int8')].';
        end
    end
end

% Repeat bits, with padding, for different CBW
switch chanBW
  case 'CBW20'
    repVHTSIGBBits = bitsPerUser; % 26
    chanBWMHz = 20;
    numCBPI   = 52;     % Table 22-6, coded bits per interleaver block
  case 'CBW40'
    repVHTSIGBBits = [bitsPerUser; bitsPerUser]; % 54
    chanBWMHz = 40;
    numCBPI   = 108;
  case 'CBW80'
    repVHTSIGBBits = [repmat(bitsPerUser, 4, 1); 0]; % 117
    chanBWMHz = 80;
    numCBPI   = 234;
  otherwise % {'CBW80+80', 'CBW160'}
    repVHTSIGBBits = repmat([repmat(bitsPerUser, 4, 1); 0], 2, 1); % 234
    chanBWMHz = 160;
    numCBPI   = 234;
end

% BCC encoding
encOut = wlanBCCEncode(repVHTSIGBBits,'1/2');

% Single symbol processing for segment deparser & interleaver
if strcmp(chanBW, 'CBW160') 
    % Two individual 80 MHz blocks
    interleaveOut = zeros(size(encOut));
    interleaveOut(1:end/2) = wlanBCCInterleave(encOut(1:2:end), ...
        'VHT', numCBPI, 1, chanBWMHz, 1);
    interleaveOut(end/2+1:end) = wlanBCCInterleave(encOut(2:2:end), ...
        'VHT', numCBPI, 1, chanBWMHz, 1);
else
    interleaveOut = wlanBCCInterleave(encOut, ...
        'VHT', numCBPI, 1, chanBWMHz, 1);
end

% Constellation mapping & segment deparser
dataPerUser = wlanConstellationMapper(interleaveOut,1);

end

% [EOF]
