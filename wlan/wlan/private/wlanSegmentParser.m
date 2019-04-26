function y = wlanSegmentParser(x, chanBW, numBPSCS, numCBPS, numES, mode)
%wlanSegmentParser Segment parser
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   Y = wlanSegmentParser(X, CHANBW, NUMBPSCS, NUMCBPS, NUMES, MODE)
%   performs segment parsing on the input X and converts it into Y when
%   CHANBW is 'CBW160'. NUMBPSCS is the number of coded bits per subcarrier
%   per spatial stream. NUMCBPS is the number of coded bits per symbol.
%   NUMES is the number of encoding streams. The MODE input can be 'Tx' or
%   'Rx', indicating the parsing is performed at transmitter or receiver,
%   respectively.
% 
%   See also wlanSegmentDeparser.

% Copyright 2015 The MathWorks, Inc.

%#codegen

% Section 22.3.10.7 in IEEE Std 802.11ac-2013
if strcmp(chanBW, 'CBW160') || strcmp(chanBW, 'CBW80+80')
    numSS     = size(x, 2); 
    numCBPSs  = numCBPS/numSS; % Input length must be multiple of numCBPSs
    blockSize = max(1, numBPSCS/2); % s in the spec
    numES     = blockSize * numES; % numES in the spec
    numTrunks = floor(numCBPSs / (2*numES));
    tailLen   = numCBPSs - 2 * numTrunks * numES;    
    
    if strcmp(mode, 'Tx')
        xIn3D = reshape(x, numCBPSs, [], numSS);  % [numCBPSs, size(x, 1)/numCBPSs, numSS]
        majorIn4D = reshape(xIn3D(1 : 2*numES*numTrunks, :, :), numES, 2, [], numSS); % [numES, 2, numTrunks*size(x,1)/numCBPSs, numSS]
        majorIn4D = permute(majorIn4D, [1 3 4 2]); % [numES, numTrunks * size(x,1)/numCBPSs, numSS, 2]

        if tailLen > 0 
            tailIn4D = reshape(xIn3D(end-tailLen+(1:tailLen), :, :), blockSize, 2, [], numSS); 
            tailIn4D = permute(tailIn4D, [1 3 4 2]);
            yIn4D = cat(1, reshape(majorIn4D, numES*numTrunks, [], numSS, 2), ...
                           reshape(tailIn4D,  tailLen/2,       [], numSS, 2)); % [obj.pnumCBPSS/2, [], numSS, 2]
            y = reshape(yIn4D, [], numSS, 2);
        else
            y = reshape(majorIn4D, [], numSS, 2);
        end    
    else
        if tailLen > 0
            xIn4D = reshape(x, numES*numTrunks+tailLen/2, [], numSS, 2);
            majorIn5D = reshape(xIn4D(1:numES*numTrunks, :, :, :), numES, numTrunks, [], numSS, 2);
            tailIn5D = reshape(xIn4D(numES*numTrunks+(1:tailLen/2), :, :, :), ...
                blockSize, (tailLen/(2*blockSize)), [], numSS, 2);
            tailIn5D = permute(tailIn5D, [1 5 2 3 4]);
        else
            majorIn5D = reshape(x, numES, numTrunks, [], numSS, 2);
            tailIn5D = [];
        end
        
        majorIn5D = permute(majorIn5D, [1 5 2 3 4]); % [numES, 2, numTrunks, size(x,1)/numCBPSs, numSS]
        if tailLen > 0
            yIn3D = cat(1, reshape(majorIn5D, 2*numES*numTrunks, [], numSS), ...
                           reshape(tailIn5D,  tailLen,           [], numSS));
        else
            yIn3D = reshape(majorIn5D, 2*numES*numTrunks, [], numSS);
        end
        
        y = reshape(yIn3D, [], numSS);
    end
else
    y = x;
end

end

