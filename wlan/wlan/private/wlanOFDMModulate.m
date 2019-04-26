function y = wlanOFDMModulate(x, CPLen)
%wlanOFDMModulate Perform OFDM modulation
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   Y = wlanOFDMModulate(X, CPLEN) performs OFDM modulation on the
%   frequency-domain input signal X and outputs the time-domain signal Y.
% 
%   X is the Nc-by-Nsym-by-Nt frequency-domain signal, where Nc represents
%   the number of frequency subcarriers, Nsym represents the number of OFDM
%   symbols, and Nt represents the number of antennas in the
%   spatial-domain.
% 
%   CPLEN is a non-negative, integer scalar representing the cyclic prefix
%   length.
% 
%   Y is the modulated ((Nc+CPLEN)*Nsym)-by-Nt time-domain signal.
%
%   See also wlanOFDMDemodulate.

% Copyright 2015 The MathWorks, Inc.

%#codegen

[FFTLen, numSym, numTx] = size(x);

% IFFT shift    
if isreal(x)
    postShift = complex(ifftshift(x, 1), 0);
else
    postShift = ifftshift(x, 1);
end

% IFFT    
postIFFT = ifft(postShift, [], 1); 
postCP = cat(1, postIFFT(end-CPLen+(1:CPLen), :, :), postIFFT);
y = reshape(postCP, [(FFTLen + CPLen)*numSym numTx]);

% [EOF]
