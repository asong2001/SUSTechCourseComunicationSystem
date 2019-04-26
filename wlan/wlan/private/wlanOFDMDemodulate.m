function [data, pilot] = wlanOFDMDemodulate(x, cfgOFDM, ofdmSymOffset)
% WLANOFDMDEMODULATE performs OFDM demodulation.
% 
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.

% Copyright 2015 The MathWorks, Inc.

%#codegen    

FFTLen    = cfgOFDM.FFTLength;
CPLen     = cfgOFDM.CyclicPrefixLength;
numRx     = size(x, 2); 
symOffset = round(ofdmSymOffset * CPLen);

% Remove cyclic prefix
if isscalar(CPLen)
    numSym = size(x, 1)/(FFTLen + CPLen);
    inputIn3D = reshape(x, [(FFTLen + CPLen) numSym numRx]);
    postCPRemoval = inputIn3D([CPLen+1:FFTLen+symOffset, symOffset+1:CPLen], :, :);
else
    numSym = length(CPLen);
    postCPRemoval = coder.nullcopy(complex(zeros(FFTLen, numSym, numRx)));
    currentIdx = 0;
    for symIdx = 1:numSym
        postCPRemoval(:, symIdx, :) = x(currentIdx + ...
            [CPLen(symIdx)+1:FFTLen+symOffset(symIdx), symOffset(symIdx)+1:CPLen(symIdx)], :);
        currentIdx = currentIdx + CPLen(symIdx) + FFTLen;
    end
end

% Denormalization
postCPRemoval = postCPRemoval / cfgOFDM.NormalizationFactor;

% FFT
postFFT = fft(postCPRemoval, [], 1);

% FFT shift
if isreal(postFFT)
    postShift = complex(fftshift(postFFT, 1), 0);
else
    postShift = fftshift(postFFT,1);
end

% Phase rotation on frequency subcarriers
postShift = bsxfun(@rdivide, postShift, cfgOFDM.CarrierRotations);

% Output data
data = postShift(cfgOFDM.DataIndices, :, :);

% Output pilots
pilot = postShift(cfgOFDM.PilotIndices, :, :);

end

% [EOF]