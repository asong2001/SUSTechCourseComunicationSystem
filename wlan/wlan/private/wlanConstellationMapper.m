function y = wlanConstellationMapper(inp, numBPSCS, varargin)
%wlanConstellationMapper Constellation mapper
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   Y = wlanConstellationMapper(INP,NUMBPSCS) maps the input bits (INP)
%   using the number of coded bits per subcarrier, NUMBPSCS. NUMBPSCS must
%   be one of 1, 2, 4, 6 or 8.
%
%   Y = wlanConstellationMapper(...,PHROT) rotates the constellation points
%   counter-clockwise by the specified amount, PHROT, in radians. This
%   applies only for NUMBPSCS = 1.
%
%   See also wlanConstellationDemodulate.

%   Copyright 2015 The MathWorks, Inc.

%#codegen
narginchk(2,3);

% Separate out BPSK from other QAM modulations
if numBPSCS==1
    % As per Fig. 18-10 Section 18.3.5.8, IEEE Std 802.11-2012
    %symMap = [0; 1];
    constellation = complex([-1; 1]);
    
    y = constellation(inp+1);
    if (nargin==3) % phRot is only used for BPSK
        phRot = varargin{1};
        y = y.*exp(1i*phRot);
    end
else
    if numBPSCS==2
        % As per Fig. 18-10 Section 18.3.5.8, IEEE Std 802.11-2012
        symMap = [1; 0; 3; 2];
    elseif numBPSCS==4
        % As per Fig. 18-10 Section 18.3.5.8, IEEE Std 802.11-2012
        symMap = [2; 3; 1; 0; 6; 7; 5; 4; 14; 15; 13; 12; 10; 11; 9; 8];
    elseif numBPSCS==6
        % As per Fig. 18-10 Section 18.3.5.8, IEEE Std 802.11-2012
        symMap = [4 5 7 6 2 3 1 0  12 13 15 14 10 11 9 8 ...
            28 29 31 30 26 27 25 24  20 21 23 22 18 19 17 16 ...
            52 53 55 54 50 51 49 48  60 61 63 62 58 59 57 56 ...
            44 45 47 46 42 43 41 40  36 37 39 38 34 35 33 32].';
    else % 256QAM
        % As per Fig. 22-24, 22-25, 22-26, 22-27, Section 22.3.10.9.1,
        % IEEE Std 802.11ac-2013
        firstCol = [8 9 11 10 14 15 13 12 4 5 7 6 2 3 1 0].';
        allQuads = [firstCol firstCol+16 firstCol+48 firstCol+32 ...
            firstCol+96 firstCol+112 firstCol+80 firstCol+64 ...
            firstCol+192 firstCol+208 firstCol+240 firstCol+224 ...
            firstCol+160 firstCol+176 firstCol+144 firstCol+128];
        symMap = allQuads(:);
    end    
        
    y = qammod(inp, 2^numBPSCS, symMap, 'InputType', 'bit', 'UnitAveragePower', true);
end

end

% [EOF]
