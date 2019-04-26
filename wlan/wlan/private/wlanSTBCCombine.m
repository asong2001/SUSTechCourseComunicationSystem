function [y, CSI] = wlanSTBCCombine(x, chanEst, numSS, eqMethod, varargin)
%wlanSTBCCombine Perform space-time block code (STBC) combining
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   [Y, CSI] = wlanSTBCCombine(X, CHANEST, NUMSS, 'ZF') performs STBC
%   combining using the signal input X and the channel estimates input
%   CHANEST, according to the STBC scheme for HT and VHT formats.
%   Zero-forcing (ZF) equalization is performed for the combining.
%
%   Y is a complex array of size Nsd-by-Nsym-by-Nss containing the estimate
%   of the transmitted signal. Nsd represents the number of data
%   subcarriers (frequency domain), Nsym represents the number of OFDM
%   symbols (time domain) and must be even, and Nss represents the number
%   of spatial streams (spatial domain), also given by NUMSS. Y is complex
%   when either X or CHANEST is complex and is real otherwise.
%
%   CSI is a real matrix of size Nsd-by-Nss representing the soft channel
%   state information.
%
%   X is complex 3D array of size Nsd-by-Nsym-by-Nr representing the
%   frequency domain input signal, where Nr represents the number of
%   receive antennas.
%
%   CHANEST is a complex 3D array of size Nsd-by-Nsts-by-Nr, where Nsts
%   represents the number of space-time streams. Nsts must be either twice
%   more than Nss or equal to Nss+1 when Nss = 2 or 3. 
%
%   NUMSS is a double precision, real, positive integer scalar, which
%   represents the number of spatial streams.
%
%   [Y, CSI] = wlanSTBCCombine(X, CHANEST, NUMSS, 'MMSE', NOISEVAR)
%   performs the STBC combining using the minimum-mean-square-error (MMSE)
%   equalization method. The noise variance input NOISEVAR is a double
%   precision, real, nonnegative scalar.
%
%   See also wlanSTBCEncode.

% Copyright 2015 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

% Input validation
narginchk(4, 5);
validateattributes(x, {'double'}, {'3d','finite','nonempty'}, ...
    'wlanSTBCCombine:InSignal', 'signal input');
coder.internal.errorIf(mod(size(x, 2), 2) ~= 0, ...
    'wlan:wlanSTBCCombine:InvalidNumSym');
validateattributes(numSS, {'double'}, {'real','scalar','positive','finite'}, ...
    'wlanSTBCCombine:NumSS', 'number of spatial streams');
validateattributes(chanEst, {'double'}, {'3d','finite','nonempty'}, ...
    'wlanSTBCCombine:ChanEst', 'channel estimation input');
coder.internal.errorIf(~strcmp(eqMethod, 'ZF') && ~strcmp(eqMethod, 'MMSE'), ...
    'wlan:wlanSTBCCombine:InvalidEqMethod');
coder.internal.errorIf(size(x, 1) ~= size(chanEst, 1), ...
    'wlan:wlanSTBCCombine:UnequalFreqCarriers');
coder.internal.errorIf(size(x, 3) ~= size(chanEst, 3), ...
    'wlan:wlanSTBCCombine:UnequalNumRx');
numSTS = size(chanEst, 2);
coder.internal.errorIf((numSTS ~= 2*numSS) && ...
    ~((numSTS == numSS + 1) && (numSS == 2 || numSS == 3)), ...
    'wlan:wlanSTBCCombine:InvalidSSAndSTSComb');

if strcmp(eqMethod, 'MMSE')
    narginchk(5,5);    
    validateattributes(varargin{1}, {'double'}, {'real','scalar','nonnegative','finite','nonempty'}, ...
        'wlanSTBCCombine:noiseVarEst', 'noise variance estimation input'); 
    noiseVarEst = varargin{1};
else
    noiseVarEst = 0;
end

% Perform combining
numSD  = size(x, 1);
numSym = size(x, 2);
numRx  = size(chanEst, 3);

% Expand channel matrix to size [numSD 2*numSS numRx] by inserting zero links
if (numSS == 2) && (numSTS == 3) 
    modifiedH = cat(2, chanEst, zeros(numSD, 1, numRx));   
elseif (numSS == 3) && (numSTS == 4)
    modifiedH = cat(2, chanEst(:,1:3,:), zeros(numSD, 1, numRx), ...
                       chanEst(:,4,:),   zeros(numSD, 1, numRx));   
else % numSTS == 2*numSS
    modifiedH = chanEst;
end
num2SS = size(modifiedH, 2); % = 2*numSS

% Derive equivalent channel matrix
eqH = complex(zeros(numSD, num2SS, 2*numRx));
eqH(:, 1:2:end, 1:2:end) =  conj(modifiedH(:, 1:2:end, :));
eqH(:, 2:2:end, 1:2:end) = -conj(modifiedH(:, 2:2:end, :));
eqH(:, :, 2:2:end) = reshape(flip(reshape(modifiedH, numSD, 2, []), 2), ...
    numSD, num2SS, numRx);

% Derive equivalent receiver matrix
x(:, 1:2:end, :) = conj(x(:, 1:2:end, :));
eqR = reshape(permute(reshape(x, numSD, 2, [], numRx), [1 3 2 4]), ...
    numSD, [], 2*numRx); % [numSD numSym/2 2*numRx]

% Perform equalization
CSI = zeros(numSD, numSS);
if (numSS == 1) && (numSTS == 2) % Take advantage of orthogonal channel matrix 
    chanEst2D = reshape(chanEst, numSD, []);
    CSI(:,1) = real(sum(chanEst2D.*conj(chanEst2D), 2)) + noiseVarEst; % [numSD numSS]
    eqOut = complex(zeros(numSD, numSym/2, numSTS)); % [numSD numSym/2 2*numSS]
    for idx = 1:numSD
        eqOut(idx, :, :) = reshape(eqR(idx,:,:), [], 2*numRx) * ...
                           reshape(eqH(idx,:,:), 2,  2*numRx)' / CSI(idx);
    end
else % Regular equalization
    % Instead of calling equalization function directly, we can exploit
    % block matrix inversion algorithm for faster performance.
    [eqOut, eqCSI] = wlanEqualize(eqR, eqH, eqMethod, noiseVarEst); % [numSD numSym/2 2*numSS]
    CSI(:) = eqCSI(:, 1:2:end); % [numSD numSS]
end

% Organize output
eqOut(:,:,1:2:end) = conj(eqOut(:,:,1:2:end));
y = reshape(permute(reshape(eqOut, numSD, [], 2, numSS), [1 3 2 4]), ...
    numSD, [], numSS); % [numSD numSym numSS]

end

% [EOF]