function out = dsssScramble(in,scramInit)
%dsssScramble DSSS scrambling
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   OUT = dsssScramble(IN,SCRAMINIT)
%

%   Copyright 2015 The MathWorks, Inc.

%#codegen

    Z = scramInit;
    
    L = length(in);
    out = zeros(L,1,'int8');
    
    for k = 1:L
        
        temp = int8(xor(Z(4),Z(7)));
        out(k) = xor(temp,in(k));
        
        Z(2:end) = Z(1:end-1);
        Z(1) = out(k);
        
    end

end

% [EOF]