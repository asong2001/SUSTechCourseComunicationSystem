function isInESSMode = inESSMode(cfgHT)
%inESSMode  Return true when NumExtensionStreams is applicable
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   ISINESSMODE = inESSMode(CFGHT) returns true when NumExtensionStreams 
%   is applicable for HT-Mixed format configuration, CFGHT.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

% cfgHT must be a wlanHTConfig object
isInESSMode = true; % by default
if ( (cfgHT.NumSpaceTimeStreams == cfgHT.NumTransmitAntennas) || ...
     cfgHT.NumSpaceTimeStreams == 4 )
    isInESSMode = false;
end

end
