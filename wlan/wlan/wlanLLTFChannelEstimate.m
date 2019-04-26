function est = wlanLLTFChannelEstimate(rxSym,cfgFormat,varargin)
% wlanLLTFChannelEstimate Channel estimation using the L-LTF
%   EST = wlanLLTFChannelEstimate(RXSYM,CFGFORMAT) returns the estimated
%   channel between the transmitter and all receive antennas using the
%   Non-HT Long Training Field (L-LTF).
%
%   EST is a complex Nst-by-1-by-Nr array containing the estimated channel
%   at data and pilot subcarriers, where Nst is the number of occupied
%   subcarriers and Nr is the number of receive antennas. The singleton
%   dimension corresponds to the single transmitted stream in the L-LTF
%   which includes the combined cyclic shifts if multiple transmit antennas
%   are used.
%
%   RXSYM is a complex Nst-by-Nsym-by-Nr array containing demodulated
%   L-LTF OFDM symbols. Nsym is the number of demodulated L-LTF symbols and
%   can be one or two. If two L-LTF symbols are provided the channel
%   estimate is averaged over the two symbols.
%
%   CFGFORMAT is a format configuration object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a>, 
%   <a href="matlab:help('wlanHTConfig')">wlanHTConfig</a>, or <a href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a>.
%
%   EST = wlanLLTFChannelEstimate(RXSYM,CHANBW) returns the estimated
%   channel using a string CHANBW, for parameterization. CHANBW is a string
%   describing the channel bandwidth which must be one of the following:
%   'CBW5','CBW10','CBW20','CBW40','CBW80','CBW160'.
%
%   EST = wlanLLTFChannelEstimate(...,SPAN) performs frequency smoothing by
%   using a moving average filter across adjacent subcarriers to reduce the
%   noise on the channel estimate. The span of the filter in subcarriers,
%   SPAN, must be odd. If adjacent subcarriers are highly correlated
%   frequency smoothing will result in significant noise reduction, however
%   in a highly frequency selective channel smoothing may degrade the
%   quality of the channel estimate. Frequency smoothing is only
%   recommended when estimating the L-LTF when a single transmit antenna is
%   used.
%
%   Examples:
%
%    Example 1:
%    %  Generate a time domain waveform txWaveform for an 802.11ac VHT
%    %  packet. Extract and demodulate the L-LTF and perform channel
%    %  estimation without frequency smoothing.
%
%       cfgVHT = wlanVHTConfig;       % Create packet configuration
%       txWaveform = wlanWaveformGenerator([1;0;0;1],cfgVHT);
%
%       rxWaveform = 0.5*txWaveform;   % Add channel gain
%       idnLLTF = wlanFieldIndices(cfgVHT,'L-LTF');
%       sym = wlanLLTFDemodulate(rxWaveform(idnLLTF(1):idnLLTF(2),:), ...
%                cfgVHT);
%       est = wlanLLTFChannelEstimate(sym,cfgVHT);
%
%    Example 2:
%    %  Generate a time domain waveform for an 802.11n HT packet, pass it
%    %  through a TGn fading channel and perform L-LTF channel estimation.
%
%       cfgHT = wlanHTConfig;         % Create packet configuration
%       txWaveform = wlanWaveformGenerator([1;0;0;1],cfgHT);
% 
%       % Configure channel
%       tgn = wlanTGnChannel;
%       tgn.SampleRate = 20e6;
% 
%       % Pass through channel (with zeros to allow for channel delay)
%       rxWaveform = step(tgn,[txWaveform; zeros(15,1)]);
%       rxWaveform = rxWaveform(5:end,:); % Synchronize for channel delay
% 
%       % Extract L-LTF and perform channel estimation
%       idnLLTF = wlanFieldIndices(cfgHT,'L-LTF');
%       sym = wlanLLTFDemodulate(rxWaveform(idnLLTF(1):idnLLTF(2),:), ...
%                cfgHT);
%       est = wlanLLTFChannelEstimate(sym,cfgHT);
%
%   See also wlanVHTConfig, wlanHTConfig, wlanNonHTConfig,
%            wlanLLTFDemodulate, wlanNonHTDataRecover, 
%            wlanHTLTFChannelEstimate, wlanVHTLTFChannelEstimate.
 
%   Copyright 2015 The MathWorks, Inc.

%#codegen

% Validate number of arguments
narginchk(2,3);

% Validate the packet format configuration object is a valid type
validateattributes(cfgFormat, ...
    {'wlanVHTConfig','wlanHTConfig','wlanNonHTConfig','char'},{}, ...
    mfilename,'first argument');

if isa(cfgFormat,'char')
    % Channel bandwidth parameterized as a char
    cbw = cfgFormat;
    coder.internal.errorIf(~(strcmpi(cbw,'CBW5') || ...
        strcmpi(cbw,'CBW10') || strcmpi(cbw,'CBW20') || ...
        strcmpi(cbw,'CBW40') || strcmpi(cbw,'CBW80')|| ...
        strcmpi(cbw,'CBW160')), ...
        'wlan:wlanChannelEstimate:InvalidChBandwidth');
else
    % Only applicable for OFDM and DUP-OFDM modulations
    coder.internal.errorIf( isa(cfgFormat, 'wlanNonHTConfig') && ...
        ~strcmp(cfgFormat.Modulation, 'OFDM'), ...
        'wlan:wlanChannelEstimate:InvalidDSSS');

    % Channel bandwidth parameterized using object
    cbw = cfgFormat.ChannelBandwidth;
end

% Validate symbol type
validateattributes(rxSym,{'single','double'},{'3d'}, ...
    'wlanLLTFChannelEstimate','L-LTF OFDM symbol(s)');
numSC = size(rxSym,1);
numRxAnts = size(rxSym,3);

% Return an empty if empty symbols
if isempty(rxSym)
    est = zeros(numSC,1,numRxAnts);
    return;
end

if nargin > 2
    span = varargin{1};
    enableFreqSmoothing = true;
else
    % Default no frequency smoothing
    enableFreqSmoothing = false;
end

% Perform LS channel estimation and time averaging as per Perahia, Eldad,
% and Robert Stacey. Next Generation Wireless LANs: 802.11 n and 802.11 ac.
% Cambridge university press, 2013, page 83, Eq 2.70.
if (strcmp(cbw,'CBW5') || strcmp(cbw,'CBW10') || strcmp(cbw,'CBW20')) 
    num20 = 1;
else
    num20 = real(str2double(cbw(4:end)))/20;
end
lltf = lltfReference(num20); % Get reference subcarriers
% Verify number of subcarriers to estimate
coder.internal.errorIf(numSC~=numel(lltf), ...
    'wlan:wlanChannelEstimate:IncorrectNumSC',numel(lltf),numSC);
ls = rxSym./repmat(lltf,1,size(rxSym,2),numRxAnts); % Least-square estimate   
est = mean(ls,2); % Average over the symbols

% Perform frequency smoothing
if enableFreqSmoothing
    % Smooth each 20 MHz segment individually
    groupSize = size(est,1)/num20;
    for i = 1:num20
        idx = (1:groupSize)+(i-1)*groupSize;
        est(idx,:,:) = frequencySmoothing(est(idx,:,:),span);
    end
end

end

function ref = lltfReference(num20MHz)
    % 20 MHz reference
    [lltfLower, lltfUpper] = lltfSequence();

    % Replicate over number of 20 MHz segments ignoring the DC and reshape
    ref = reshape([lltfLower; lltfUpper]*ones(1,num20MHz),[],1);
end
