function y = wlanSegmentDeparser(x, chanBW, mode)
%wlanSegmentDeparser Segment deparser
%
%   Note: This is an internal undocumented function and its
%   API and/or functionality may change in subsequent releases.
%
%   Y = wlanSegmentDeparser(X, CHANBW, MODE) performs segment deparsing on
%   the input X and converts it into Y when CHANBW is 'CBW160'. The MODE
%   input can be 'Tx' or 'Rx', indicating the deparsing is performed at
%   transmitter or receiver, respectively.
% 
%   See also wlanSegmentParser.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

% Section 22.3.10.9.3 in IEEE Std 802.11ac-2013
if strcmp(chanBW, 'CBW160') 
    NSD80 = 234;
    if strcmp(mode, 'Tx')
        numSS = size(x, 2);
        tempX = reshape(x, NSD80, [], numSS, 2); 
        tempX = permute(tempX, [1 4 2 3]); 
        y = reshape(tempX, [], numSS);             
    else % x is of dimension [Nsd, Nsym, numSS]
        Nsym = size(x, 2);
        numSS  = size(x, 3);
        tempX = reshape(x, NSD80, 2, [], Nsym*numSS);
        tempX = permute(tempX, [1 3 4 2]);
        y = reshape(tempX, [], Nsym, numSS, 2);        
    end
else
    y = x;
end

end

