function y = wlanVHTSTF(cfgVHT)
%WLANVHTSTF VHT Short Training Field (VHT-STF)
% 
%   Y = wlanVHTSTF(CFGVHT) generates the VHT Short Training Field (VHT-STF)
%   time-domain signal for the VHT transmission format.
%
%   Y is the time-domain VHT-STF signal. It is a complex matrix of size
%   Ns-by-Nt where Ns represents the number of time-domain samples and Nt
%   represents the number of transmit antennas.
%
%   CFGVHT is the format configuration object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a> which
%   specifies the parameters for the VHT format.
% 
%   Example:
%   %  Generate the VHT-STF signal for a VHT 80MHz transmission format
% 
%     cfgVHT = wlanVHTConfig;                % Format configuration
%     cfgVHT.ChannelBandwidth = 'CBW80';     % Set to 80MHz 
%     vstfOut = wlanVHTSTF(cfgVHT);
% 
%   See also wlanVHTConfig, wlanLSTF, wlanVHTLTF.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

validateattributes(cfgVHT, {'wlanVHTConfig'}, {'scalar'}, mfilename, ...
                   'VHT format configuration object');
validateConfig(cfgVHT, 'SMapping');
               
chanBW = cfgVHT.ChannelBandwidth;
numSTSTotal = sum(cfgVHT.NumSpaceTimeStreams);    

% Get OFDM parameters
cfgOFDM = wlanGetOFDMConfig(chanBW, 'Long', 'VHT', numSTSTotal);
FFTLen = cfgOFDM.FFTLength;
CPLen  = cfgOFDM.CyclicPrefixLength;
num20  = FFTLen/64;

% Non-HT L-STF (IEEE Std:802.11-2012, pg 1695)
VHTSTF = lstfSequence();

numVHTFtones = 12*num20;  % Defined as per Table 22-8 (page 252)
vhtf = [zeros(6,1);  VHTSTF; zeros(5,1)];

% Replicate over CBW and apply phase rotation
vhtfToneRotated = repmat(vhtf, num20, 1) .* cfgOFDM.CarrierRotations;

% Replicate over multiple antennas
vhtfMIMO = repmat(vhtfToneRotated, 1, numSTSTotal);

% Cyclic shift addition.
% The cyclic shift is applied per space-time stream
csh = getCyclicShiftVal('VHT', numSTSTotal, 20*num20);
vhtfCycShift = wlanCyclicShift(vhtfMIMO, csh, FFTLen, 'Tx');

% Spatial mapping
vhtfSpatialMapped = wlanSpatialMapping(vhtfCycShift, cfgVHT.SpatialMapping, ...
    cfgVHT.NumTransmitAntennas, cfgVHT.SpatialMappingMatrix);

% OFDM modulation
modOut = ifft(ifftshift(vhtfSpatialMapped, 1), [], 1);

out = [modOut; modOut(1:CPLen,:)];  % Extend to a long symbol duration
y = out.*(FFTLen/sqrt(numVHTFtones*numSTSTotal));

end

% [EOF]
