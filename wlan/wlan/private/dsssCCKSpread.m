function waveform = dsssCCKSpread(phi)
%dsssCCKSpread DSSS CCK spreading
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   WAVEFORM = dsssCCKSpread(PHI)
%

%   Copyright 2015 The MathWorks, Inc.

%#codegen

    % Each row of 'argc' defines the combinations of phase terms 
    % phi1...phi4 for each chip of the CCK codeword C
    % Clause 17.4.6.6.1 General
    argc = [ 1  1  1  1;
             1  0  1  1;
             1  1  0  1;
             1  0  0  1;
             1  1  1  0;
             1  0  1  0;
             1  1  0  0;
             1  0  0  0; ];

    % Each row of 'magc' defines the magnitude for each chip of the CCK
    % codeword C
    % Clause 17.4.6.6.1 General
    % Clause 17.4.6.6.2 Cover code for CCK
    magc = [1 1 1 -1 1 1 -1 1].';
    
    % Establish cckSymbols, the number of CCK symbols in the input
    cckSymbols = size(phi,1);

    % Columns of 'c' are the overall CCK codewords
    c = repmat(magc,1,cckSymbols).*exp(1i*argc*phi.');

    % Reshape 'c' into the overall waveform
    waveform = c(:);

end

% [EOF]