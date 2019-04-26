function startIdx = helperSymbolTiming(x, chanBW, varargin)
%helperSymbolTiming Perform symbol timing based on L-LTF in time-domain
%   STARTIDX = helperSymbolTiming(X, CHANBW) returns the estimated starting
%   index for the Non-HT Long Training Field (L-LTF), using
%   cross-correlation method.
%   
%   STARTIDX is an integer scalar indicating the starting index relative to
%   X for the first L-LTF this function successfully detects. STARTIDX
%   could be a non-positive number, which typically happens when X starts
%   in the middle of a L-LTF. STARTIDX is set to empty when this function
%   does not detect any L-LTF.
%   
%   X is the received time-domain signal on which symbol timing is
%   performed. It is a Ns-by-Nr matrix of real or complex values, where Ns
%   represents the number of time-domain samples and Nr represents the
%   number of receive antennas.
% 
%   CHANBW is the channel bandwidth, one of 'CBW20', 'CBW40', 'CBW80' and
%   'CBW160'.
%   
%   STARTIDX = helperSymbolTiming(..., THRESHOLD) specifies the relative
%   threshold for L-LTF peak detection after cross-correlation as a real
%   scalar between 0 and 1, exclusive. When unspecified a value of 0.3 is
%   used.
%
%   % Example: 
%   %    Introduce a delay to L-LTF and perform symbol timing on it. 
%   
%   cfgHT  = wlanHTConfig('ChannelBandwidth', 'CBW40');
%   txLLTF = wlanLLTF(cfgHT);           % L-LTF generation
%   txSig  = [zeros(15,1); txLLTF];     % Introduce a delay
%   rxSig  = awgn(txSig, 1, 1);         % Add noise
%   
%   % Perform symbol timing and expect the returned index to be 16
%   LLTFIdx = helperSymbolTiming(rxSig, 'CBW40')
%
%   % If the input starts in L-LTF, the returned index can be 0 or negative
%   LLTFIdx = helperSymbolTiming(rxSig(17:end,:), 'CBW40') % Should be 0
%   LLTFIdx = helperSymbolTiming(rxSig(18:end,:), 'CBW40') % Should be -1 
% 
%   See also wlanCoarseCFOEstimate, wlanFineCFOEstimate.

%#codegen

narginchk(2,3);

validateattributes(x, {'double'}, {'2d','finite'}, mfilename, 'signal input'); 
if isempty(x)
    startIdx = [];
    return;
end

if ~any(strcmp(chanBW,{'CBW5','CBW10','CBW20','CBW40','CBW80','CBW160'}))
    error(['The chanBW input must be one of ''CBW5'', ''CBW10'',' ...
           ' ''CBW20'', ''CBW40'', ''CBW80'' and ''CBW160''.']);
end

if nargin == 2
    threshold = .3;
else
    validateattributes(varargin{1}, {'double'}, ...
        {'real','scalar','>',0,'<',1}, mfilename, 'threshold');
    threshold = varargin{1};
end

FFTLen = helperFFTLength(chanBW);
symLen = FFTLen * 5/4;
bufLen = 4*symLen;
xLen = size(x, 1);

if (strcmp(chanBW,'CBW5') || strcmp(chanBW,'CBW10')) 
    chanBW = 'CBW20'; % override for VHTConfig use below.
end
cfgVHT = wlanVHTConfig('ChannelBandwidth', chanBW);
LLTF = wlanLLTF(cfgVHT);
coeff = conj(flipud(LLTF(1:symLen,1)));

if xLen <= FFTLen
    startIdx = [];
elseif xLen <= bufLen
    corr = sum(abs(filter(coeff, 1, x)).^2, 2);
    peakIdx = find(corr > max(corr)*threshold);
    startIdx = processCorr(peakIdx, FFTLen, symLen, bufLen, 0);
else % xLen > bufLen
    bufEndIdx = 0;
    filterState = complex(zeros(symLen-1, size(x, 2)));
    numSampNext = bufLen;
    corr = zeros(bufLen, 1);
    startIdxInBuf = [];
    while numSampNext > 0
        [corrNew, filterState] = filter(coeff, 1, x(bufEndIdx+(1:numSampNext), :), filterState);
        bufEndIdx = bufEndIdx + numSampNext;
        corr(1:bufLen) = [corr(numSampNext+1:bufLen, 1); sum(abs(corrNew).^2, 2)];
        peakIdx = find(corr > max(corr)*threshold);
        [startIdxInBuf, numSampNext] = processCorr(peakIdx, FFTLen, symLen, bufLen, xLen - bufEndIdx);
    end    
    startIdx = bufEndIdx - bufLen + startIdxInBuf;
end

startIdx = startIdx - symLen + 1;

end

function [startIdx, numSampNext] = processCorr(peakIdx, FFTLen, symLen, bufLen, leftSamples)

numPeaks = length(peakIdx);

if numPeaks > .1*bufLen % Most likely noise input if 10% samples are over threshold
    startIdx = [];
    numSampNext = min(symLen, leftSamples);
else % do some validation
    [numPairs, lastPair1stIdx, lastPair2ndIdx] = deal(0);
    i = 1;  secondClusterStartIdx = numPeaks+1;
    while i < secondClusterStartIdx 
        thisPair2ndIdx = find(peakIdx(i) + FFTLen == peakIdx);
        if ~isempty(thisPair2ndIdx)
            numPairs = numPairs + 1;
            lastPair1stIdx = peakIdx(i);
            lastPair2ndIdx = peakIdx(thisPair2ndIdx(1));
            if (numPairs == 1) && (secondClusterStartIdx == numPeaks + 1)
                secondClusterStartIdx = thisPair2ndIdx(1);
            end        
        end
        i = i + 1;      
    end
    
    if numPairs == 0 
        startIdx = [];
        numSampNext = min(3*symLen, leftSamples); 
    elseif (lastPair2ndIdx > bufLen - symLen/2) && (leftSamples > 0) % numPairs > 0
        startIdx = [];
        numSampNext = min(symLen, leftSamples); 
    else
        startIdx = lastPair1stIdx;
        numSampNext = 0; 
    end
end

end

% [EOF]