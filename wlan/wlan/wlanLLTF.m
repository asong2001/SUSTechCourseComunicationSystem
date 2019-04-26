function y = wlanLLTF(cfgFormat)
%WLANLLTF Non-HT Long Training Field (L-LTF)
%
%   Y = wlanLLTF(CFGFORMAT) generates the Non-HT Long Training Field
%   (L-LTF) time-domain signal for the VHT, HT-Mixed, and Non-HT OFDM
%   transmission formats.
%
%   Y is the time-domain L-LTF signal. It is a complex matrix of size
%   Ns-by-Nt, where Ns represents the number of time-domain samples and
%   Nt represents the number of transmit antennas.
%
%   CFGFORMAT is the format configuration object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a>, 
%   <a href="matlab:help('wlanHTConfig')">wlanHTConfig</a>, or <a href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a>, which specifies the parameters for
%   the VHT, HT-Mixed, and Non-HT OFDM formats, respectively. Only OFDM 
%   modulation is supported for a wlanNonHTConfig object input.
%
%   Example: Generate the L-LTF signal for a VHT 80MHz SISO packet
%
%     cfgVHT = wlanVHTConfig('ChannelBandwidth', 'CBW80');
%     y = wlanLLTF(cfgVHT);
%     plot(abs(y));
%
%   See also wlanVHTConfig, wlanHTConfig, wlanNonHTConfig, wlanLSTF, wlanLSIG.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

% Validate the format configuration object
validateattributes(cfgFormat, {'wlanVHTConfig','wlanHTConfig','wlanNonHTConfig'}, ...
    {'scalar'}, mfilename, 'format configuration object');

% Only applicable for OFDM and DUP-OFDM modulations
coder.internal.errorIf( isa(cfgFormat, 'wlanNonHTConfig')  && ...
                        ~strcmp(cfgFormat.Modulation, 'OFDM'), ...
                        'wlan:wlanLLTF:InvalidNonHTLLTF');

chanBW = cfgFormat.ChannelBandwidth;
if (strcmp(chanBW, 'CBW10') || strcmp(chanBW, 'CBW5'))
    numTx = 1;  % override cfgFormat and set to 1 only, for 802.11j/p
else
    numTx = cfgFormat.NumTransmitAntennas;
end
 
% OFDM params
cfgOFDM = wlanGetOFDMConfig(chanBW, 'Half', 'Legacy');
FFTLen = cfgOFDM.FFTLength;
num20  = FFTLen/64;

[lltfLower, lltfUpper] = lltfSequence();
lltf = [zeros(6,1);  lltfLower; 0; lltfUpper; zeros(5,1)];

% Replicate over channel bandwidth & Tx, and apply phase rotation
lltfMIMO = bsxfun(@times, repmat(lltf, num20, numTx), cfgOFDM.CarrierRotations);

% Cyclic shift addition.
% The FORMAT is set to OFDM due to legacy mode. The cyclic shift is
% applied on each transmit antenna.
csh = getCyclicShiftVal('OFDM', numTx, 20*num20);
lltfCycShift = wlanCyclicShift(lltfMIMO, csh, FFTLen, 'Tx');

modOut = ifft(ifftshift(lltfCycShift, 1), [], 1);        

% CP length = TGI2  (IEEE Std 802.11ac-2013, Table 22-8)
out = [modOut((end-FFTLen/2+1):end,:); modOut; modOut];

y  = out * cfgOFDM.NormalizationFactor / sqrt(numTx);

end

% [EOF]
