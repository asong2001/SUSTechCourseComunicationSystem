function y = wlanBCCDecode(x, rate, varargin)
%wlanBCCDecode Binary convolutional coding decoder.
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   Y = wlanBCCDecode(X,RATE) convolutionally decodes the input
%   (X) using a WLAN convolutional code at the specified rate (RATE).
%   RATE must be one of:
%   '5/6' or 5/6, which maps to a puncture pattern of [1 1 1 0 0 1 1 0 0 1]
%   '3/4' or 3/4, which maps to a puncture pattern of [1 1 1 0 0 1]
%   '2/3' or 2/3, which maps to a puncture pattern of [1 1 1 0]
%   '1/2' or 1/2, which maps to no puncturing
%
%   Y = wlanBCCDecode(...,TBD) allows for the traceback depth to be
%   specified as TBD.
%
%   See also wlanBCCEncode.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

% Add soft-decision capability - hardcode to 3 bits.
% Validate the 3 inputs
%   Rate - enum with 4 strings or double
%   x - input vector
%   tbd - scalar integer

narginchk(2,3);
assert(isreal(x));
if nargin==3
    tbd = varargin{1};
    applyTR = false;     % applyThumbRule
else
    applyTR = true;
    tbd = 30; % dummy
end

trellis = poly2trellis(7, [133 171]);
puncPat = [1 1]; % to be overridden below

% Traceback depth thumb rule = 2.5m/(1-r), where m=6
% Override the thumb rule:
%   for SIGB decoding in 20MHZ = 26.
%   for L-SIG decoding = 24.
if ~strcmp(rate, '1/2') && ~isequal(rate, 1/2)
    if strcmp(rate, '2/3') || isequal(rate, 2/3)
        puncPat = [1 1 1 0].';
        if applyTR, tbd = 45; end
        y = zeros(round(size(x,1)*2/3), 1, 'int8');
    elseif strcmp(rate, '3/4') || isequal(rate, 3/4)
        puncPat = [1 1 1 0 0 1].';
        if applyTR, tbd = 60; end
        y = zeros(round(size(x,1)*3/4), 1, 'int8');
    else  % strcmp(rate, '5/6') || isequal(rate, 5/6)
        puncPat = [1 1 1 0 0 1 1 0 0 1].';
        if applyTR, tbd = 90; end
        y = zeros(round(size(x,1)*5/6), 1, 'int8');
    end
    
    y(:) = int8(vitdec(x, trellis, tbd, 'trunc', 'unquant', puncPat));
else % rate = '1/2'
    if applyTR, tbd = 30; end    
    y = zeros(round(size(x,1)/2), 1, 'int8');
    y(:) = int8(vitdec(x, trellis, tbd, 'trunc', 'unquant'));
end

end

% [EOF]
