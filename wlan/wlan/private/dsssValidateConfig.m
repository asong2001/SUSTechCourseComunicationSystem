function dsssValidateConfig(cfgNonHT,caller)
%dsssValidateConfig validation of DSSS configuration
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   dsssValidateConfig(CFGNONHT,CALLER)
%

%   Copyright 2015 The MathWorks, Inc.

%#codegen

    validateattributes(cfgNonHT,{'wlanNonHTConfig'},{'scalar'},mfilename,'format configuration object');
    coder.internal.errorIf(~strcmp(cfgNonHT.Modulation,'DSSS'),['wlan:' caller ':InvalidModulation']);
    
end

% [EOF]