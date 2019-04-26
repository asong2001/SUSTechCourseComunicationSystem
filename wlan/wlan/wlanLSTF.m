function y = wlanLSTF(cfgFormat)
%WLANLSTF Non-HT Short Training Field (L-STF)
%
%   Y = wlanLSTF(CFGFORMAT) generates the Non-HT Short Training Field
%   (L-STF) time-domain signal for the VHT, HT-Mixed, and Non-HT OFDM
%   transmission formats.
%
%   Y is the time-domain L-STF signal. It is a complex matrix of size
%   Ns-by-Nt, where Ns represents the number of time-domain samples and
%   Nt represents the number of transmit antennas.
%
%   CFGFORMAT is the format configuration object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a>, 
%   <a href="matlab:help('wlanHTConfig')">wlanHTConfig</a>, or <a href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a>, which specifies the parameters for
%   the VHT, HT-Mixed, and Non-HT OFDM formats, respectively. Only OFDM 
%   modulation is supported for a wlanNonHTConfig object input.
%
%   Example: Generate L-STF for a single transmit antenna VHT 80MHz packet
%
%     cfgVHT = wlanVHTConfig('ChannelBandwidth', 'CBW80');
%     y = wlanLSTF(cfgVHT);
%     plot(abs(y));
%
%   See also wlanVHTConfig, wlanHTConfig, wlanNonHTConfig, wlanLLTF.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

% Validate the format configuration object
validateattributes(cfgFormat, {'wlanVHTConfig','wlanHTConfig','wlanNonHTConfig'}, ...
    {'scalar'}, mfilename, 'format configuration object');

% Only applicable for OFDM and DUP-OFDM modulations
coder.internal.errorIf( isa(cfgFormat, 'wlanNonHTConfig') && ...
                        ~strcmp(cfgFormat.Modulation, 'OFDM'),...
                        'wlan:wlanLSTF:InvalidNonHTLSTF');

chanBW = cfgFormat.ChannelBandwidth;
if (strcmp(chanBW, 'CBW10') || strcmp(chanBW, 'CBW5'))
    numTx = 1;  % override cfgFormat and set to 1 only, for 802.11j/p
else
    numTx = cfgFormat.NumTransmitAntennas;
end

% OFDM parameters
cfgOFDM = wlanGetOFDMConfig(chanBW, 'Half', 'Legacy');
FFTLen = cfgOFDM.FFTLength;
num20  = FFTLen/64;

% Non-HT L-STF (IEEE Std:802.11-2012, pg 1695)
LSTF = lstfSequence();

% The short training field consists of 12 subcarriers
N_LSTF_TONE = 12*num20;
lstf = [zeros(6,1); LSTF; zeros(5,1)];

% Replicate over channel bandwidth & Tx, and apply phase rotation
lstfMIMO = bsxfun(@times, repmat(lstf, num20, numTx), cfgOFDM.CarrierRotations);

% Cyclic shift addition.
% The FORMAT is set to OFDM due to legacy mode. The cyclic shift is
% applied on each transmit antenna.
csh = getCyclicShiftVal('OFDM', numTx, 20*num20);
lstfCycShift = wlanCyclicShift(lstfMIMO, csh, FFTLen, 'Tx');

modOut = ifft(ifftshift(lstfCycShift, 1), [], 1);       
out = [modOut; modOut; modOut(1:FFTLen/2,:)];

y = out * (FFTLen/sqrt(numTx*N_LSTF_TONE));

end

% [EOF]
