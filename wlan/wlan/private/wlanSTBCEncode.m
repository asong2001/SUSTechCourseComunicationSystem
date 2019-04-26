function y = wlanSTBCEncode(x, numSTS)
%wlanSTBCEncode Perform space-time block coding (STBC)
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   Y = wlanSTBCEncode(X, NUMSTS) encodes the input X according to the STBC
%   scheme for HT and VHT formats, and returns the result in Y. The input X
%   can be a double precision 2-D matrix or 3-D array with real or complex
%   values. X is of size Nsd x Nsym x Nss, where Nsd represents the number
%   of data subcarriers (frequency domain), Nsym represents the number of
%   OFDM symbols (time domain) that must be even, and Nss represents the
%   number of spatial streams (spatial domain). The input NUMSTS is a
%   double precision, real, positive integer scalar, which represents the
%   number of space-time streams Nsts. Nsts must be either twice more than
%   Nss or equal to Nss+1 when Nss = 2 or 3. The output Y is of size Nsd x
%   Nsym x Nsts. Y has the same data type and complexity as X.
%
%   See also wlanSTBCCombine.

% Copyright 2015 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

% Input validation
validateattributes(x, {'double'}, {'3d','finite','nonempty'}, ...
    'wlanSTBCEncode:InSignal', 'signal input');
coder.internal.errorIf(mod(size(x, 2), 2) ~= 0, ...
    'wlan:wlanSTBCEncode:InvalidNumSym');
validateattributes(numSTS, {'double'}, {'real','scalar','finite','>=',2}, ...
    'wlanSTBCEncode:NumSTS', 'number of space-time streams');
numSS = size(x, 3);
coder.internal.errorIf((numSTS ~= 2*numSS) && ...
    ~((numSTS == numSS + 1) && (numSS == 2 || numSS == 3)), ...
    'wlan:wlanSTBCEncode:InvalidSSAndSTSComb');

% Perform encoding: Section 20.3.11.9.2 in IEEE Std 802.11-2012 and Section
% 22.3.10.9.4 in IEEE Std 802.11ac-2013
if (numSS == 2 && numSTS == 3) || (numSS == 3 && numSTS == 4) % Specific HT
    stsStayIdx = [1:2:numSS, numSTS];
    ssFlipIdx  = 1;
else % numSTS == 2*numSS for VHT and HT
    stsStayIdx = 1:2:numSTS;
    ssFlipIdx  = 1:numSS;
end
stsFlipIdx = 2*ssFlipIdx;

y = complex(zeros(size(x, 1), size(x, 2), numSTS));
y(:,:,stsStayIdx) = x;
y(:,:,stsFlipIdx) = reshape(bsxfun(@times, fliplr(...
    conj(reshape(x(:,:,ssFlipIdx), size(x, 1), 2, []))), [-1 1]), size(x(:,:,ssFlipIdx)));

% [EOF]