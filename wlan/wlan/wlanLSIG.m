function [y, bits] = wlanLSIG(cfgFormat)
%WLANLSIG Non-HT SIGNAL field (L-SIG)
%
%   [Y, BITS] = wlanLSIG(CFGFORMAT) generates the Non-HT SIGNAL field
%   (L-SIG) time-domain signal for the VHT, HT-Mixed, and Non-HT OFDM
%   transmission formats.
%
%   Y is the time-domain L-SIG field signal. It is a complex matrix of size
%   Ns-by-Nt, where Ns represents the number of time-domain samples and
%   Nt represents the number of transmit antennas.
%
%   BITS is the signaling bits used for the Non-HT SIGNAL field. It is
%   an int8-typed, binary column vector of length 24.
%
%   CFGFORMAT is the format configuration object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a>, 
%   <a href="matlab:help('wlanHTConfig')">wlanHTConfig</a>, or <a href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a>, which specifies the parameters for
%   the VHT, HT-Mixed, and Non-HT OFDM formats, respectively. Only OFDM 
%   modulation is supported for a wlanNonHTConfig object input.
%
%   Example: 
%   % Generate the LSIG signal for a VHT 80MHz transmission format
%
%     cfgVHT = wlanVHTConfig;              % Format configuration
%     cfgVHT.ChannelBandwidth = 'CBW80';   % Set to 80MHz
%     y = wlanLSIG(cfgVHT);
%
%   See also wlanVHTConfig, wlanHTConfig, wlanNonHTConfig, wlanLLTF, 
%   wlanLSIGRecover.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

% Validate the format configuration object
validateattributes(cfgFormat, {'wlanVHTConfig','wlanHTConfig', ...
    'wlanNonHTConfig'}, {'scalar'}, 'wlanLSIG', ...
    'format configuration object');

if isa(cfgFormat, 'wlanVHTConfig')

    % Validate VHT configuration
    s = validateConfig(cfgFormat, 'MCSSTSTx');
    txTime = s.TxTime;
    
    % Set the RATE value to 4 bit binary code. The code value is fixed
    % to [1 1 0 1] representing 6Mb/s in legacy 20MHz CBW.
    % As per 22.3.8.2.4, p/256, IEEE Std 802.11ac-2013.
    R = [1; 1; 0; 1];
    length = ceil((txTime - 20)/4)*3 - 3;

elseif isa(cfgFormat, 'wlanHTConfig')

    % Validate HT configuration
    s = validateConfig(cfgFormat, 'MCSSTSTx');
    txTime = s.TxTime;
        
    % Set the RATE value to 4 bit binary code. The code value is fixed
    % to [1 1 0 1] representing 6Mb/s in legacy 20MHz CBW.
    % As per Sec. 20.2.9.3.6 and Sec. 9.23.4, IEEE Std 802.11-2012
    R = [1; 1; 0; 1];
    % Assuming signalExtension = 0, from Sec 9.23.4.
    length = ceil((txTime - 20)/4)*3 - 3;

elseif isa(cfgFormat, 'wlanNonHTConfig')
    % Only applicable for OFDM and DUP-OFDM modulations
    coder.internal.errorIf( ~strcmp(cfgFormat.Modulation, 'OFDM'), ...
                            'wlan:wlanLSIG:InvalidNonHTLSIG');
    validateConfig(cfgFormat);

    switch cfgFormat.MCS
        case 0 % 6 Mbps
            R = [1; 1; 0; 1];
        case 1 % 9 Mbps
            R = [1; 1; 1; 1];
        case 2 % 12 Mbps
            R = [0; 1; 0; 1];
        case 3 % 18 Mbps
            R = [0; 1; 1; 1];
        case 4 % 24 Mbps
            R = [1; 0; 0; 1];
        case 5 % 36 Mbps
            R = [1; 0; 1; 1];
        case 6  % 48 Mbps
            R = [0; 0; 0; 1];
        otherwise % 7 => 54 Mbps
            R = [0; 0; 1; 1];
    end
    
    length = cfgFormat.PSDULength;

end
chanBW = cfgFormat.ChannelBandwidth;

% Construct the SIGNAL field
% Length parameter with LSB first, which is 12 bits 
lengthBits  = de2bi(length, 12, 'right-msb').';

% Even parity bit 
parityBit = mod(sum([R;lengthBits],1),2);

% The SIGNAL field (As per IEEE Std 802.12-2012, 18.3.4, pg. 1595)
bits = [R; 0; lengthBits; parityBit; zeros(6, 1, 'int8')];

%% Process L-SIG bits
encodedBits = wlanBCCEncode(bits, '1/2');
interleavedBits = wlanBCCInterleave(encodedBits, 'NON_HT', 48, 1);
modData = wlanConstellationMapper(interleavedBits, 1);

% OFDM prms based on bandwidth
% Use long cyclic prefic length as per Table 22-8, IEEE Std 802.11ac-2013
cfgOFDM = wlanGetOFDMConfig(chanBW, 'Long', 'Legacy');
FFTLen = cfgOFDM.FFTLength;
CPLen  = cfgOFDM.CyclicPrefixLength;
num20  = FFTLen/64;

% Add pilot symbols, from IEEE Std 802.11-2012, Eqn 20-14
dataSymbol = complex(zeros(FFTLen, 1));
dataSymbol(cfgOFDM.DataIndices, 1) = repmat(modData, num20, 1);
Nsym = 1; % One symbol
z = 0;    % No offset as first symbol with pilots
dataSymbol(cfgOFDM.PilotIndices, 1) = repmat(nonHTPilots(Nsym, z), num20, 1);

% Phase rotation
lsigToneRotated = dataSymbol .* cfgOFDM.CarrierRotations;

% Replicate over multiple antennas
if (strcmp(chanBW, 'CBW10') || strcmp(chanBW, 'CBW5'))
    numTx = 1;  % override and set to 1 only, for 802.11j/p
else
    numTx = cfgFormat.NumTransmitAntennas;
end
lsigMIMO = repmat(lsigToneRotated(:), 1, numTx);

% Cyclic shift addition.
% The FORMAT is set to OFDM due to legacy mode. The cyclic shift is
% applied on each transmit antenna.
csh = getCyclicShiftVal('OFDM', numTx, 20*num20);
lsigCycShift = wlanCyclicShift(lsigMIMO, csh, FFTLen, 'Tx');

modOut = ifft(ifftshift(lsigCycShift, 1), [], 1);
out = [modOut((end-CPLen+1):end,:); modOut]; % Add cyclic prefix

y = out * cfgOFDM.NormalizationFactor / sqrt(numTx);

end

% [EOF]
