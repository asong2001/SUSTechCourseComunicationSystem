function estSmooth = frequencySmoothing(est,span)
%frequencySmoothing Moving average filtering across subcarriers

%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   ESTSMOOTH = frequencySmoothing(EST,SPAN) perform moving average
%   filtering in the frequency domain across subcarriers using an odd
%   length window. The window span is decreased at the ends as to not bias
%   the result. 
%
%   ESTSMOOTH is a Ns-by-Nsts-by-Nr matrix containing the smoothed channel
%   estimates, where Ns is the number of subcarriers, Nsts is the number of
%   space-time streams and Nr is the number of receive antennas.
%
%   EST is a Ns-by-Nsts-by-Nr matrix containing the channel estimates to be
%   smoothed.
%
%   SPAN is the smoothing span to use. 

%   Copyright 2015 The MathWorks, Inc.

%#codegen

% Validate smoothing span attributes
validateattributes(span,{'numeric'},{'>=',1,'odd','scalar',} ...
    ,'','smoothingSpan');

% No smoothing when the span is 1
if span==1
    estSmooth = est;
    return;
end
span = cast(span,'double');

numSC = size(est,1);
numSTS = size(est,2);
numRxAnts = size(est,3);

estSmooth = complex(zeros(size(est)),zeros(size(est)));
for sts = 1:numSTS
    for rx = 1:numRxAnts
        weights = ones(span,1); % Use a rectangular window
        data = [est(:,sts,rx); zeros(span-1,1)];
        averageData = filter(weights,1,data);

        % Remove unwanted elements
        removeIdx = [2:2:span (length(data)-(2:2:span)+1) ...
            numSC:2:((span-1)*(mod(numSC,2)==1))]; % For span>numSym
        averageData(removeIdx) = [];

        % Normalization factor given that the effective window size
        % changes at the edges when not enough subcarriers are
        % available
        M = min(numSC-1,span-2);
        MM = min(numSC-~mod(numSC,2),span-2);
        normFactor = [1:2:M span*ones(1,length(averageData)-(span-1)) MM:-2:1].';

        estSmooth(:,sts,rx) = averageData./normFactor;
    end
end
end