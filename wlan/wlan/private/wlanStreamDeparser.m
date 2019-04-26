function y = wlanStreamDeparser(x, Nes, Nbpscs)  
%wlanStreamDeparser stream deparser.
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%

% Copyright 2015 The MathWorks, Inc.

%#codegen 

if (size(x, 2) == 1) && (Nes == 1)
    y = x;
else
    Nss        = size(x, 2);
    blockSize  = max(1, Nbpscs/2);
    inLen      = size(x, 1);
    tailLen    = rem(inLen, blockSize * Nes);
    
    tempX = reshape(x(1:inLen - tailLen, :), blockSize, Nes, [], Nss);
    tempX = permute(tempX, [1 4 3 2]);
    y     = reshape(tempX, [], Nes);

    if tailLen > 0
        tempX2 = reshape(x(end-tailLen+(1:tailLen), :), blockSize, [], Nss);
        tempX2 = permute(tempX2, [1 3 2]); 
        y = [y; reshape(tempX2, [], Nes)];
    end    
end

end

% [EOF]