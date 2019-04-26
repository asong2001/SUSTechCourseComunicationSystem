function [y, err] = wlanCRCDetect(x)
%wlanCRCDetect CRC detection.
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%

% Copyright 2015 The MathWorks, Inc.

%#codegen

y = x(1:end-8, :);
checksum = wlanCRCGenerate(y);
err = any(checksum ~= x(end-8+(1:8), :));
