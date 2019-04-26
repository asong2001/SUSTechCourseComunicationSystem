function y = wlanSpatialMapping(x, mappingType, numTx, mappingMatrix)
%wlanSpatialMapping Spatial mapping
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   Y = wlanSpatialMapping(X, MAPPINGTYPE, NUMTX, MAPPINGMATRIX) performs
%   spatial mapping from space-time streams to transmit antennas. X is a
%   FFTLen-by-numSTS or numST-by-numSTS matrix, where FFTLen represents the
%   FFT length, numST represents the number of data plus pilot subcarriers,
%   and numSTS represents the number of space-time streams. MAPPINGTYPE can
%   be one of 'Direct', 'Hadamard', 'Fourier' and 'Custom'. NUMTX is the
%   number of transmit antennas. MAPPINGMATRIX is a numSTS-by-NUMTX,
%   FFTLen-by-numSTS-by-NUMTX or numST-by-numSTS-by-NUMTX spatial mapping
%   matrix(ces) that apply only when the MAPPINGTYPE input is 'Custom'. The
%   output Y is a FFTLen-by-NUMTX or numST-by-NUMTX matrix.

% Copyright 2015 The MathWorks, Inc.

%#codegen

numCarriers = size(x, 1); % = FFTLen or numST
numSTS      = size(x, 2); % = numSTS

% Section 20.3.10.11.1 in IEEE Std 802.11-2012.
switch mappingType
  case 'Direct'
    y = x;
  case 'Hadamard'
    dataAndPilotIdx = getDataAndPilotIndices(numCarriers);
    y = complex(zeros(numCarriers, numTx));
    Q = hadamard(8);
    normQ = Q(1:numSTS, 1:numTx)/sqrt(numTx);
    y(dataAndPilotIdx, :) = x(dataAndPilotIdx, :) * normQ;
  case 'Fourier'
    dataAndPilotIdx = getDataAndPilotIndices(numCarriers);
    y = complex(zeros(numCarriers, numTx));
    % The following can be obtained from dftmtx(numTx) which however does not generate code
    [g1, g2] = meshgrid(0:numTx-1, 0:numSTS-1);
    normQ = exp(-1i*2*pi.*g1.*g2/numTx)/sqrt(numTx);
    y(dataAndPilotIdx, :) = x(dataAndPilotIdx, :) * normQ;
  otherwise  % case 'Custom'
    dataAndPilotIdx = getDataAndPilotIndices(numCarriers);
    y = complex(zeros(numCarriers, numTx));
    if size(mappingMatrix, 1) <= 8
        Q = mappingMatrix(1:numSTS, :);
        normQ = Q * sqrt(numSTS)/norm(Q, 'fro'); % Normalization
        y(dataAndPilotIdx, :) = x(dataAndPilotIdx, :) ...
            * normQ(:, 1:numTx); % Need the 1:numTx indexing for codegen
    else
        for idx = 1:length(dataAndPilotIdx)
            freqIdx = dataAndPilotIdx(idx);            
            Q = reshape(mappingMatrix(idx, 1:numSTS, :), numSTS, numTx);
            normQ = Q * sqrt(numSTS)/norm(Q, 'fro');            
            y(freqIdx, :) = x(freqIdx, :) * normQ;
        end
    end
end

end

function dataAndPilotIdx = getDataAndPilotIndices(numCarriers)

if any(numCarriers == [56 114 242 484])
    dataAndPilotIdx = 1:numCarriers;
else
    FFTLen = numCarriers;
    DCOffset = FFTLen/2 + 1;
    switch FFTLen
      case 64
        nullIndices = [1:4 DCOffset FFTLen-2:FFTLen];
      case 512
        nullIndices = [1:6 DCOffset + [-129:-127 -5:5 127:129] FFTLen-4:FFTLen]; 
      otherwise
        nullIndices = [1:6 DCOffset + (-1:1) FFTLen-4:FFTLen];
    end

    dataAndPilotIdx = setdiff(1:FFTLen, nullIndices);
end

end