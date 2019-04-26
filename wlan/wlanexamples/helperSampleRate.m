function fs = helperSampleRate(cfg)
% helperSampleRate Return the nominal sample rate for format configuration
%
%   FS = helperSampleRate(CFGFORMAT) returns the nominal sample rate for
%   the specified format configuration object, CFGFORMAT.
%
%   FS is a scalar representing the sample rate in Hz.
%
%   CFGFORMAT is the format configuration object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a>, 
%   <a href="matlab:help('wlanHTConfig')">wlanHTConfig</a>, or <a href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a>, which specifies the parameters for
%   the VHT, HT-Mixed, and Non-HT formats, respectively.
%
%   FS = helperSampleRate(CHANBW) returns the nominal sample rate for the
%   specified channel bandwidth. CHANBW must be one of 'CBW5', 'CBW10',
%   'CBW20', 'CBW40', 'CBW80' or 'CBW160'.
%
%   Example: Return sample rate for a VHT format configuration. 
%
%   cfgVHT = wlanVHTConfig;
%   fs = helperSampleRate(cfgVHT)
%
%   fs = helperSampleRate(cfgVHT.ChannelBandwidth)
%
%   See also wlanVHTConfig, wlanHTConfig, wlanNonHTConfig.

%#codegen

if ischar(cfg)
    coder.internal.errorIf( ~any(strcmp(cfg, {'CBW5', 'CBW10', ...
            'CBW20', 'CBW40', 'CBW80', 'CBW160'})), ...
        'wlan:helperSampleRate:InvalidChBandwidth');
    
    fs = real(str2double(cfg(4:end)))*1e6;    
else    
    validateattributes(cfg, {'wlanVHTConfig','wlanHTConfig', ...
        'wlanNonHTConfig'}, {'scalar'}, mfilename, ...
        'format configuration object');

    if ( isa(cfg, 'wlanVHTConfig') || isa(cfg, 'wlanHTConfig') || ...
        (isa(cfg, 'wlanNonHTConfig') && strcmp(cfg.Modulation, 'OFDM')) )
        fs = real(str2double(cfg.ChannelBandwidth(4:end)))*1e6;

    else % NonHT DSSS
        fs = 11e6;
    end
end

% [EOF]
