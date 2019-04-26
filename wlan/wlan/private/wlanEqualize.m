function [y, CSI] = wlanEqualize(x, chanEst, eqMethod, varargin)
%wlanEqualize Perform MIMO channel equalization. 
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   [Y, CSI] = wlanEqualize(X, CHANEST, 'ZF') performs equalization using
%   the signal input X and the channel estimation input CHANEST, and
%   returns the estimation of transmitted signal in Y and the soft channel
%   state information in CSI. The zero-forcing (ZF) method is used. The
%   inputs X and CHANEST can be double precision 2-D matrices or 3-D arrays
%   with real or complex values. X is of size Nsd x Nsym x Nr, where Nsd
%   represents the number of data subcarriers (frequency domain), Nsym
%   represents the number of OFDM symbols (time domain), and Nr represents
%   the number of receive antennas (spatial domain). CHANEST is of size Nsd
%   x Nsts x Nr, where Nsts represents the number of space-time streams.
%   The double precision output Y is of size Nsd x Nsym x Nsts. Y is
%   complex when either X or CHANEST is complex and is real otherwise. The
%   double precision, real output CSI is of size Nsd x Nsts.
%
%   [Y, CSI] = wlanEqualize(X, CHANEST, 'MMSE', NOISEVAR) performs the
%   equalization using the minimum-mean-square-error (MMSE) method. The
%   noise variance input NOISEVAR is a double precision, real, nonnegative
%   scalar.
%
%   See also wlanSTBCCombine.

% Copyright 2015 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

% Input validation
narginchk(3, 4);

validateattributes(x, {'double'}, {'3d','finite','nonempty'}, ...
    'wlanEqualize:InSignal', 'signal input');
validateattributes(chanEst, {'double'}, {'3d','finite','nonempty'}, ...
    'wlanEqualize:ChanEst', 'channel estimation input');   
coder.internal.errorIf(~strcmp(eqMethod, 'ZF') && ~strcmp(eqMethod, 'MMSE'), ...
    'wlan:wlanEqualize:InvalidEqMethod');
coder.internal.errorIf(size(x, 1) ~= size(chanEst, 1), ...
    'wlan:wlanEqualize:UnequalFreqCarriers');
coder.internal.errorIf(size(x, 3) ~= size(chanEst, 3), ...
    'wlan:wlanEqualize:UnequalNumRx');

if strcmp(eqMethod, 'MMSE')
    narginchk(4,4);
    validateattributes(varargin{1}, {'double'}, {'real','scalar','nonnegative','finite','nonempty'}, ...
        'wlanEqualizer:noiseVarEst', 'noise variance estimation input'); 
    noiseVarEst = varargin{1};
else % ZF
    noiseVarEst = 0;
end

% Perform equalization
numTx  = size(chanEst, 2);
numRx  = size(chanEst, 3);

CSI = zeros(size(x, 1), numTx); % Pre-allocation here for code generation
if (numTx == 1 && numRx == 1) % SISO
    CSI = real(chanEst.*conj(chanEst)) + noiseVarEst;
    y =  bsxfun(@times, x, conj(chanEst(:))./CSI(:)); % (:) for codegen
elseif (numTx == 1 && numRx > 1) % SIMO
    chanEst2D = reshape(chanEst, size(chanEst, 1), numRx);    
    CSI = real(diag(chanEst2D*chanEst2D')) + noiseVarEst;
    y = bsxfun(@rdivide, sum(bsxfun(@times, x, conj(chanEst)), 3), CSI);
elseif (numTx > 1 && numRx == 1) % MISO
    chanPower = real(chanEst.*conj(chanEst));
    if strcmp(eqMethod, 'ZF')
        CSI = chanPower; 
    else % Use Schur complement formula
        CSI = noiseVarEst + noiseVarEst*chanPower./...
            (bsxfun(@minus, sum(chanPower, 2), chanPower) + noiseVarEst);
    end
    chanEstInv = bsxfun(@rdivide, conj(chanEst), sum(chanPower, 2)+noiseVarEst);
    y = bsxfun(@times, x(:,:,1), permute(chanEstInv, [1 3 2])); % (:,:,1) for codegen
elseif (numTx > numRx) && strcmp(eqMethod, 'ZF') % MIMO: singular channel matrix using ZF
    numSym = size(x, 2);
    CSI = sum(real(chanEst .* conj(chanEst)), 3);
    y = complex(zeros(size(x, 1), numSym, numTx));
    for idx = 1:size(chanEst, 1)
        y(idx, :, 1:numTx) = reshape(x(idx, :, :), numSym, numRx) * ...
                  pinv(reshape(chanEst(idx,:,:), numTx, numRx)); 
    end
else % MIMO: numTx > numRx using MMSE or numTx <= numRx
    numSym = size(x, 2);
    y = complex(zeros(size(x, 1), numSym, numTx));
    for idx = 1:size(chanEst, 1)
        H = reshape((chanEst(idx,:,:)), numTx, numRx);
        invH = inv(H*H'+noiseVarEst*eye(numTx));
        CSI(idx, :)  = 1./real(diag(invH));
        y(idx, :, 1:numTx) = reshape(x(idx, :, :), numSym, numRx) * H' * invH;  %#ok<MINV>
    end    
end

% [EOF