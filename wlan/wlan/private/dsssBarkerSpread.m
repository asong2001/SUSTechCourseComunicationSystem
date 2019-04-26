function waveform = dsssBarkerSpread(symbols)
%dsssBarkerSpread DSSS barker spreading
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   WAVEFORM = dsssBarkerSpread(SYMBOLS)
%

%   Copyright 2015 The MathWorks, Inc.

%#codegen

    barker = [1 -1  1  1 -1  1  1  1 -1 -1 -1].';
    
    waveform = kron(symbols,barker);

end

% [EOF]