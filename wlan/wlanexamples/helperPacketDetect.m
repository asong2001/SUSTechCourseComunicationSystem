function startIndex = helperPacketDetect(x, chanBW, varargin)
%HELPERPACKETDETECT detects the start of a received OFDM WLAN packet
%
%   STARTINDEX = helperPacketDetect(X,CHANBW) returns the start index of
%   the preamble of the received WLAN packet, using auto-correlation
%   method.
%
%   STARTINDEX is an integer scalar indicating the start index of a
%   detected WLAN packet relative to X. If the packet is not detected an
%   empty value is returned indicating an undetected packet.
%
%   X is a received real or complex time-domain signal of size Ns-by-Nr
%   matrix, where Ns is the number of time domain samples and Nr is the
%   number of receive antennas.
% 
%   CHANBW is a string describing the channel bandwidth and must be 'CBW5',
%   'CBW10', 'CBW20' 'CBW40', 'CBW80' or 'CBW160'.
%
%   STARTIDX = helperPacketDetect(..., THRESHOLD) specifies the absolute
%   threshold for L-STF peak detection after auto-correlation as a real
%   scalar greater than 0 and less than or equal to 1. When unspecified a
%   value of 0.5 is used.
%
%   Example:
%   %  Detect a received 802.11n packet at 20dB SNR.
%
%      cfgHT = wlanHTConfig();    % Create packet configuration
%      SNR = 20;                  % In decibels
%      tgn = wlanTGnChannel('LargeScaleFadingEffect','None');
%   
%      % Generate transmit waveform
%      txWaveform = wlanWaveformGenerator([1;0;0;1],cfgHT);
% 
%      fadedSig = step(tgn,txWaveform); 
%      rxWaveform = awgn(fadedSig, SNR, 0);
%
%      startIndex = helperPacketDetect(rxWaveform,cfgHT.ChannelBandwidth);
%
%   See also helperSymbolTiming, wlanFieldIndices.
 
%   Copyright 2015 The MathWorks, Inc.

%#codegen

narginchk(2,3);
validateattributes(x, {'double'}, {'2d','finite'}, mfilename, 'signal input'); 

if isempty(x)
    startIndex = [];
    return;
end

if nargin == 2
    threshold = 0.5;
else
    validateattributes(varargin{1}, {'double'}, ...
        {'real','scalar','>',0,'<=',1}, mfilename, 'threshold');
    threshold = varargin{1};
end

% Validate channel bandwidth
coder.internal.errorIf(~any(strcmp(chanBW,{'CBW5','CBW10','CBW20', ...
    'CBW40','CBW80','CBW160'})), ...
    'wlan:helperPacketDetect:InvalidChBandwidth');

Td = 0.8e-6; % Time period of a short training symbol for 20MHz
    
switch chanBW
    case 'CBW40'
        symbolLength = Td/(1/40e6);
    case {'CBW80'}
        symbolLength = Td/(1/80e6);
    case {'CBW160'}
        symbolLength = Td/(1/160e6);
    otherwise % 'CBW5', 'CBW10', 'CBW20'
        symbolLength = Td/(1/20e6); % In samples
end
lstfLength = symbolLength*10; % Length of 10 L-STF symbols 

coder.internal.errorIf(size(x,1) < lstfLength, ...
            'wlan:helperPacketDetect:InvalidInputLength');

inpLength = size(x,1);
% Pad input with zeros if the size is less than L-STF length
if mod(inpLength,lstfLength)
    padSamples = zeros(lstfLength-mod(inpLength,lstfLength),size(x,2));
    x = [x;padSamples];
end

% Process the input waveform in blocks of L-STF length. The processing
% blocks are offset by half the L-STF length.
numBlocks = inpLength/(lstfLength/2)-1;
startIndex = [];

for n=0:numBlocks-1
    % Update buffer
    buffer = x(n*lstfLength/2+(1:lstfLength),:);
    startIndex = correlateSamples(buffer,symbolLength,lstfLength,threshold);
    
    if ~isempty(startIndex)
        % Packet detected
        startIndex = startIndex + n*lstfLength/2;
        break;
    end

end

end

function packetStart = correlateSamples(rxSig,symbolLength,lstfLength,threshold)
%   Estimate the start index of the preamble of the receive WLAN packet,
%   using auto-correlation method [1,2].
%
%   [1] OFDM Wireless LANs: A Theoretical and Practical Guide 1st Edition
%       by Juha Heiskala (Author),John Terry Ph.D. ISBN-13:978-0672321573
%   [2] OFDM Baseband Receiver Design for Wireless Communications by
%       Tzi-Dar Chiueh, Pei-Yun Tsai. ISBN: 978-0-470-82234-0

    correlationLength = lstfLength -(symbolLength*2)+1;
    pNoise = eps; % Addiing noise to avoid the divide by zero
    weights = ones(symbolLength,1);
    index = 1:correlationLength;

    % Initialize output
    packetStart = [];

    % Shift data for correlation
    rxDelayed = rxSig(symbolLength+1:end,:); % Old samples
    rx = rxSig(1:end-symbolLength,:);        % New samples
    
    % Sum output on multiple receive antennas
    C = sum(filter(weights,1,(conj(rxDelayed).*rx)),2);
    CS = C(symbolLength:end)./symbolLength;
    
    % Sum output on multiple receive antennas
    P = sum(filter(weights,1,(abs(rxDelayed).^2 + ...
                        abs(rx).^2)/2)./symbolLength,2);
    PS = P(symbolLength:end) + pNoise;

    M = abs(CS).^2 ./ PS.^2;

    N = M > threshold;
 
    if (sum(N) >= symbolLength*1.5)
        found = index(N);
        packetStart = found(1);  % Pick first peak
        
        % Check the relative distance between peaks relative to the first
        % peak. If this exceed three times the symbol length then the
        % packet is not detected.
        if sum((found(2:symbolLength) - found(1))> symbolLength*3)
            packetStart = [];
        end
    end
end

% [EOF]
