function y = wlanStreamParser(x, numSS, numBPSCS)  
%wlanStreamParser Stream parser
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   Y = wlanStreamParser(X, NUMSS, NUMBPSCS) performs stream parsing on the
%   input X and converts it into Y. NUMSS is the number of spatial streams.
%   NUMBPSCS is the number of coded bits per subcarrier per spatial stream.
%   X and Y have Nes and NUMSS columns respectively, where Nes represents
%   the number of encoding streams.
%
%   See also wlanStreamDeparser.

% Copyright 2015 The MathWorks, Inc.

%#codegen 

if (size(x, 2) == 1) && (numSS == 1)
    y = x;
else
    % Section 20.3.11.8.2 in in IEEE Std 802.11-2012 and Section
    % 22.3.10.6 in IEEE Std 802.11ac-2013.
    numES   = size(x, 2);
    blkSize = max(1, numBPSCS/2);
    inLen   = size(x, 1);
    tailLen = rem(inLen, blkSize*numSS); % One segment is equal to blkLen * numSS
    
    tempX = reshape(x(1:inLen - tailLen, :), blkSize, numSS, [], numES); % [blkSize, numSS, numSegs, numES]
    tempX = permute(tempX, [1 4 3 2 ]); % [blkSize, numES, numSegs, numSS]
    y = reshape(tempX, [], numSS);
    
    if tailLen > 0
        tempX2 = reshape(x(end-tailLen+(1:tailLen), :), blkSize, numSS, []);
        tempX2 = permute(tempX2, [1 3 2]);
        y = [y; reshape(tempX2, [], numSS)];
    end
end

end

