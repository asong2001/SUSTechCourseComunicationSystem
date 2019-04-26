function y = wlanLLTFDemodulate(rxLLTF,cfgFormat,varargin)
%WLANLLTFDEMODULATE OFDM demodulate L-LTF signal
%
%   Y = wlanLLTFDemodulate(RXLLTF,CHANBW) demodulates the time-domain
%   Non-HT Long training field (L-LTF) received signal for VHT, HT-Mixed,
%   and Non-HT OFDM transmission formats.
%
%   Y is the frequency-domain signal corresponding to the L-LTF.
%   It is a complex matrix or 3-D array of size Nst-by-2-by-Nr, where Nst
%   represents the number of used subcarriers in the L-LTF, and Nr
%   represents the number of receive antennas. Two OFDM symbols are
%   demodulated for the L-LTF.
%
%   RXLLTF is the received time-domain L-LTF signal. It is a complex matrix
%   of size Ns-by-Nr, where Ns represents the number of time-domain
%   samples. Ns can be greater than or equal to the L-LTF length, lenLLTF,
%   where only the first lenLLTF samples of RXLLTF are used.
%
%   CHANBW is a string representing the channel bandwidth, one of 'CBW5',
%   'CBW10', 'CBW20', 'CBW40', 'CBW80', or 'CBW160'.
%
%   Y = wlanLLTFDemodulate(RXLLTF,CFGFORMAT) similarly demodulates the
%   received signal when a transmission configuration object is specified.
%   CFGFORMAT is of type <a href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a> or <a href="matlab:help('wlanHTConfig')">wlanHTConfig</a> or <a
%   href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a> 
%   that specifies the channel bandwidth. Only OFDM modulation is supported
%   for a wlanNonHTConfig object input.
%
%   Y = wlanLLTFDemodulate(...,SYMOFFSET) specifies the optional OFDM
%   symbol sampling offset as a fraction of the cyclic prefix length
%   between 0 and 1, inclusive. When unspecified a value of 0.75 is used.
%
%   Example: Demodulate a received L-LTF signal 
%
%     cfgVHT = wlanVHTConfig;                   % VHT format configuration
%     txLLTF = wlanLLTF(cfgVHT);                % L-LTF generation
%       
%     rxLLTF = awgn(txLLTF, 1, 1);              % Add noise
%     y = wlanLLTFDemodulate(rxLLTF,cfgVHT);    % Demodulate
%
%   See also wlanLLTF, wlanVHTConfig, wlanLLTFChannelEstimate. 

%   Copyright 2015 The MathWorks, Inc.

%#codegen

narginchk(2,3);

% cfgFormat validation - string or object
if ischar(cfgFormat) % BW input
    coder.internal.errorIf( ~( strcmpi(cfgFormat,'CBW5') || ...
        strcmpi(cfgFormat,'CBW10') || ...
        strcmpi(cfgFormat,'CBW20') || ...
        strcmpi(cfgFormat,'CBW40') || ...
        strcmpi(cfgFormat,'CBW80') || ...
        strcmpi(cfgFormat,'CBW160')), ...
        'wlan:wlanLLTFDemodulate:InvalidChBandwidth');
    chanBW = cfgFormat;
else
    % Validate the format configuration object
    validateattributes(cfgFormat, {'wlanVHTConfig','wlanHTConfig', ...
        'wlanNonHTConfig'}, {'scalar'}, mfilename, ...
        'format configuration object');
    
    % Only applicable for OFDM and DUP-OFDM modulations
    coder.internal.errorIf( isa(cfgFormat, 'wlanNonHTConfig') && ...
        ~strcmp(cfgFormat.Modulation, 'OFDM'), ...
        'wlan:wlanLLTFDemodulate:InvalidDSSS');

    chanBW = cfgFormat.ChannelBandwidth;        
end

% Input rxLLTF validation
validateattributes(rxLLTF, {'double'}, {'2d', 'finite'}, ...
    'rxLLTF', 'L-LTF signal'); 

if nargin == 3
    validateattributes(varargin{1}, {'double'}, ...
        {'real','scalar','>=',0,'<=',1}, mfilename, 'symOffset');
    symOffset = varargin{1};
else    % default
    symOffset = 0.75;
end

numRx = size(rxLLTF, 2);
if size(rxLLTF, 1) == 0
    y = zeros(0, 0, numRx);
    return;
end

cfgOFDM = wlanGetOFDMConfig(chanBW, 'Half', 'Legacy');
FFTLen = cfgOFDM.FFTLength;
[~, sortedDataPilotIdx] = sort([cfgOFDM.DataIndices; cfgOFDM.PilotIndices]);

% Validate length of input
coder.internal.errorIf(size(rxLLTF, 1) < 2.5*FFTLen, ...
    'wlan:wlanLLTFDemodulate:ShortDataInput', 2.5*FFTLen);

% OFDM demodulation
[ofdmDemodData, ofdmDemodPilots] = wlanOFDMDemodulate([rxLLTF(1:1.5*FFTLen, :); ...
                        rxLLTF(FFTLen+(1:1.5*FFTLen), :)], cfgOFDM, symOffset);
                    
% Sort data and pilot subcarriers
ofdmDemod = [ofdmDemodData; ofdmDemodPilots];
y = ofdmDemod(sortedDataPilotIdx, :, :);

end
