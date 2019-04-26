function y = wlanConstellationDemodulate(inp,numBPSCS,nVar,varargin)
%wlanConstellationDemodulate Constellation demodulation.
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   Y = wlanConstellationDemodulate(INP, NUMBPSCS, NVAR) demodulates the
%   received input symbols (INP) using the soft-decision approximate LLR
%   method for the specified number of coded bits per subcarrier, NUMBPSCS.
%   NUMBPSCS must be one of 1, 2, 4, 6 or 8.. INP must be a 2D or 3D double
%   precision array where the last dimension specifies the number of
%   received streams. NVAR specifies the noise variance per received stream
%   and can be a scalar or vector.
%
%   Y = wlanConstellationDemodulate(..., PHROT) derotates the symbols by
%   the specified amount in radians, PHROT. This applies only for NUMBPSCS
%   = 1.
%
%   See also wlanConstellationMapper.

%   Copyright 2015 The MathWorks, Inc.

%#codegen
narginchk(3,4);

% Validate the input symbols
validateattributes(inp, {'double'}, {'nonempty', 'nonsparse', ...
    'finite'}, mfilename, 'input symbols'); 

% Validate nVar wrt to input
if ~isscalar(nVar) % is vector, check length
    inpSize = size(inp);
    % nVar applies along last dimension of input
    coder.internal.errorIf(length(nVar) ~= inpSize(end), ...
    'wlan:wlanConstellationDemodulate:InvalidNVarInput');
end

% Clip nVar to allowable value to avoid divide by zero warnings
minNVar = 1e-10;
if any(nVar < minNVar)
    nVar(nVar < minNVar) = minNVar;
end

% Separate out BPSK from other QAM modulations
if numBPSCS==1   
    if (nargin==4) % phRot is only used for BPSK, derotate
        phRot = varargin{1};
        inp = inp.*exp(-1i*phRot);
    end    
    
    % For constellation [-1, 1] as per Fig. 18-10 Section 18.3.5.8, 802.11-2012
    approxLLR = -4*real(inp); % = abs(1 - inp).^2 - abs(-1 - inp).^2;

    if isscalar(nVar) % scalar n0 applies to entire input
        y = approxLLR / nVar;
    else
        inpSize = size(inp); 
        nVarSize = [ones(1, length(inpSize)-1), inpSize(end)];
        y = bsxfun(@rdivide, approxLLR, reshape(nVar, nVarSize));
    end
else
    if numBPSCS==2
        % As per Fig. 18-10 Section 18.3.5.8, 802.11-2012
        symMap = [1; 0; 3; 2];
    elseif numBPSCS==4
        % As per Fig. 18-10 Section 18.3.5.8, 802.11-2012
        symMap = [2; 3; 1; 0; 6; 7; 5; 4; 14; 15; 13; 12; 10; 11; 9; 8];
    elseif numBPSCS==6
        % As per Fig. 18-10 Section 18.3.5.8, 802.11-2012
        symMap = [4 5 7 6 2 3 1 0  12 13 15 14 10 11 9 8 ...
            28 29 31 30 26 27 25 24  20 21 23 22 18 19 17 16 ...
            52 53 55 54 50 51 49 48  60 61 63 62 58 59 57 56 ...
            44 45 47 46 42 43 41 40  36 37 39 38 34 35 33 32].';
    else % 256QAM
        % As per Fig. 22-24, 22-25, 22-26, 22-27, Section 22.3.10.9.1,
        % 11ac Std.
        firstCol = [8 9 11 10 14 15 13 12 4 5 7 6 2 3 1 0].';
        allQuads = [firstCol firstCol+16 firstCol+48 firstCol+32 ...
            firstCol+96 firstCol+112 firstCol+80 firstCol+64 ...
            firstCol+192 firstCol+208 firstCol+240 firstCol+224 ...
            firstCol+160 firstCol+176 firstCol+144 firstCol+128];
        symMap = allQuads(:);
    end
    
    y = qamdemod(inp, 2^numBPSCS, symMap, 'UnitAveragePower', true, ...
        'OutputType', 'approxllr', 'NoiseVariance', nVar);        
end

end

% [EOF]
