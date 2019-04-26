function y = wlanHTSTF(cfgHT)
%wlanHTSTF HT Short Training Field (HT-STF)
%
%   Y = wlanHTSTF(CFGHT) generates the HT Short Training Field (HT-STF)
%   time-domain signal for the HT-Mixed transmission format.
%
%   Y is the time-domain HT-STF signal. It is a complex matrix of size
%   Ns-by-Nt, where Ns represents the number of time-domain samples and
%   Nt represents the number of transmit antennas.
%
%   CFGHT is the format configuration object of type <a href="matlab:help('wlanHTConfig')">wlanHTConfig</a> which
%   specifies the parameters for the HT-Mixed format.
%
%   Example: 
%   %  Generate the HT-STF signal for a HT 40MHz transmission format
%
%      cfgHT = wlanHTConfig;                % Format configuration
%      cfgHT.ChannelBandwidth = 'CBW40';    % Set to 40MHz
%      htstfOut = wlanHTSTF(cfgHT);
%
%   See also wlanHTConfig, wlanLSTF, wlanHTSIG, wlanHTLTF.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

validateattributes(cfgHT, {'wlanHTConfig'}, {'scalar'}, mfilename, ...
                   'HT-Mixed format configuration object');
validateConfig(cfgHT, 'SMapping'); 

numSTS = cfgHT.NumSpaceTimeStreams;

% OFDM parameters
cfgOFDM = wlanGetOFDMConfig(cfgHT.ChannelBandwidth, 'Long', 'HT', numSTS);
FFTLen = cfgOFDM.FFTLength;
CPLen  = cfgOFDM.CyclicPrefixLength;
num20  = FFTLen/64;

% Non-HT L-STF (IEEE Std:802.11-2012, pg 1695)
HTSTF = lstfSequence();

NHTFtones = 12*num20;  % as per Table 20-8, Std 802.11-2012
htf = [zeros(6,1);  HTSTF; zeros(5,1)];

% Replicate over channel bandwidth & numSTS, and apply phase rotation
htfMIMO = bsxfun(@times, repmat(htf, num20, numSTS), cfgOFDM.CarrierRotations);

% Cyclic shift applied per STS
csh = getCyclicShiftVal('VHT', numSTS, 20*num20);
htfCycShift = wlanCyclicShift(htfMIMO, csh, FFTLen, 'Tx');

% Spatial mapping
htfSpatialMapped = wlanSpatialMapping(htfCycShift, cfgHT.SpatialMapping, ...
    cfgHT.NumTransmitAntennas, cfgHT.SpatialMappingMatrix);

% OFDM modulation
modOut = ifft(ifftshift(htfSpatialMapped, 1), [], 1);

out = [modOut; modOut(1:CPLen,:)];
y = out * (FFTLen/sqrt(NHTFtones*numSTS));

end

% [EOF]
