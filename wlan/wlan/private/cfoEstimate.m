function foffset = cfoEstimate(x,D)
%cfoEstimate Frequency offset estimation
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   FOFFSET = cfoEstimate(X,D) estimates the carrier frequency offset
%   FOFFSET in Hertz using a time-domain sequence containing a repeated
%   pattern.
%
%   X is a complex Ns-by-Nr matrix containing the time domain signal. Ns is
%   the number of time domain samples, and Nr is the number of receive
%   antennas. The signal contains a repeated pattern with a period D. At
%   least two repetitions of the sequence must be present.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

% At least two repetitions of the pattern are required for estimation
[numSamples,numRxAnts] = size(x);
coder.internal.errorIf(numSamples<(2*D), ...
'wlan:wlanCFOEstimate:NotEnoughSamples');

% CFO estimate with multiple receive antennas
% Van Zelst and Schenk, Implementation of a MIMO OFDM-Based
% Wireless LAN System, IEEE TRANSACTIONS ON SIGNAL PROCESSING, VOL.
% 52, NO. 2, FEBRUARY 2004
unused = mod(numSamples,D); % Unused samples at the end
cx = x(1:end-(D+unused),:);
sx = x(D+1:(end-unused),:);

% Calculate for each receive antenna
res = complex(zeros(numRxAnts,1));
for n = 1:numRxAnts
    res(n) = cx(:,n)'*sx(:,n);
end
foffset = angle(sum(res))/(2*pi);

end