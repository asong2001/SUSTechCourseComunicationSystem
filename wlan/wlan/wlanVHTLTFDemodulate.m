function y = wlanVHTLTFDemodulate(rxVHTLTF,cfgVHT,varargin)
%wlanVHTLTFDemodulate OFDM demodulate VHT-LTF signal
%
%   Y = wlanVHTLTFDemodulate(RXVHTLTF,CFGVHT) demodulates the time-domain
%   VHT-LTF received signal for the VHT transmission format.
%   
%   Y is the frequency-domain signal corresponding to the VHT-LTF. It is a
%   complex matrix or 3-D array of size Nst-by-Nsym-by-Nr, where Nst
%   represents the number of data and pilot subcarriers in the VHT-LTF,
%   Nsym represents the number of OFDM symbols in the VHT-LTF, and Nr
%   represents the number of receive antennas.
%
%   RXVHTLTF is the received time-domain VHT-LTF signal. It is a complex
%   matrix of size Ns-by-Nr, where Ns represents the number of samples. Ns
%   can be greater than or equal to the VHT-LTF length, lenVHT, where only
%   the first lenVHT samples of RXVHTLTF are used.
%   
%   CFGVHT is the format configuration object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a>, which
%   specifies the parameters for the VHT format.
%
%   Y = wlanVHTLTFDemodulate(...,SYMOFFSET) specifies the optional OFDM
%   symbol sampling offset as a fraction of the cyclic prefix length
%   between 0 and 1, inclusive. When unspecified, a value of 0.75 is used.
%
%   Example: Demodulate a received VHT-LTF signal 
%
%     cfgVHT = wlanVHTConfig;                   % VHT format configuration
%     txVHTLTF = wlanVHTLTF(cfgVHT);            % VHT-LTF generation
%       
%     rxVHTLTF = awgn(txVHTLTF, 1, 1);            % Add noise
%     y = wlanVHTLTFDemodulate(rxVHTLTF,cfgVHT);  % Demodulate
%
%   See also wlanVHTLTF, wlanVHTConfig, wlanVHTLTFChannelEstimate.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

narginchk(2,3);

% cfgVHT validation
validateattributes(cfgVHT, {'wlanVHTConfig'}, {'scalar'}, mfilename, ...
                   'VHT format configuration object');
% Dependent validation not needed for necessary fields (CHANBW, numSTS)
coder.internal.errorIf(cfgVHT.NumUsers ~= 1, ...
    'wlan:wlanVHTLTFDemodulate:SingleUserOnly');

% Input rxVHTLTF validation
validateattributes(rxVHTLTF, {'double'}, {'2d', 'finite'}, ...
    'rxVHTLTF', 'VHT-LTF signal'); 

numRx = size(rxVHTLTF, 2);
if size(rxVHTLTF, 1) == 0
    y = zeros(0, 0, numRx);
    return;
end

if nargin == 3
    validateattributes(varargin{1}, {'double'}, ...
        {'real','scalar','>=',0,'<=',1}, mfilename, 'SYMOFFSET');
    symOffset = varargin{1};
else    % default
    symOffset = 0.75;
end

numSymTable = [1, 2, 4, 4, 6, 6, 8, 8];
numSTSTotal = sum(cfgVHT.NumSpaceTimeStreams);
numSym = numSymTable(numSTSTotal);

% Get OFDM configuration
cfgOFDM = wlanGetOFDMConfig(cfgVHT.ChannelBandwidth, 'Long', ...
    'VHT', numSTSTotal);
[~, sortedDataPilotIdx] = sort([cfgOFDM.DataIndices; cfgOFDM.PilotIndices]);
    
% Cross-validation between inputs
minInpLen = numSym*(cfgOFDM.FFTLength+cfgOFDM.CyclicPrefixLength);
coder.internal.errorIf(size(rxVHTLTF, 1) < minInpLen, ...
    'wlan:wlanVHTLTFDemodulate:ShortDataInput', minInpLen);
    
% OFDM demodulation
[ofdmDemodData, ofdmDemodPilots] = ...
    wlanOFDMDemodulate(rxVHTLTF(1:minInpLen,:), cfgOFDM, symOffset);

% Sort data and pilot subcarriers
ofdmDemod = [ofdmDemodData; ofdmDemodPilots];
y = ofdmDemod(sortedDataPilotIdx, :, :);

end

% [EOF]
