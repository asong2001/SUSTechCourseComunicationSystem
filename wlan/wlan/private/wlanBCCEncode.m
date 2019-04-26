function y = wlanBCCEncode(x,rate)
%wlanBCCEncode Binary convolutional coding encoder.
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   Y = wlanBCCEncode(X, RATE) convolutionally encodes the binary input
%   data X using a WLAN convolutional code at the specified rate (RATE).
%   RATE must be one of:
%   '5/6' or 5/6, which maps to a puncture pattern of [1 1 1 0 0 1 1 0 0 1]
%   '3/4' or 3/4, which maps to a puncture pattern of [1 1 1 0 0 1]
%   '2/3' or 2/3, which maps to a puncture pattern of [1 1 1 0]
%   '1/2' or 1/2, which maps to no puncturing
%
%   See also wlanBCCDecode.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

% Expected inputs
%   Rate - enum with 4 strings or double
%   x - input binary vector

% Reference: IEEE Std 802.11-2012, Sections 18.3.5.6, 18.3.5.7, 20.3.11.6 
% for polynomials and puncturing patterns.

trellis = poly2trellis(7, [133 171]);
if ~strcmp(rate, '1/2') && ~isequal(rate, 1/2)
    if strcmp(rate, '2/3') || isequal(rate, 2/3)
        puncPat = [1; 1; 1; 0];
        y = zeros(round(size(x,1)*3/2), 1, 'int8');
    elseif strcmp(rate, '3/4') || isequal(rate, 3/4)
        puncPat = [1; 1; 1; 0; 0; 1];
        y = zeros(round(size(x,1)*4/3), 1, 'int8');
    else % (Rate, '5/6')
        puncPat = [1; 1; 1; 0; 0; 1; 1; 0; 0; 1];
        y = zeros(round(size(x,1)*6/5), 1, 'int8');
    end
    
    y(:) = int8(convenc(x, trellis, puncPat));
else % rate 1/2 - no puncturing    
    y = zeros(size(x,1)*2, 1, 'int8');
    y(:) = int8(convenc(x, trellis));    
end

end

% [EOF]
