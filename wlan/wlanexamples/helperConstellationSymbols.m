function const = helperConstellationSymbols(cfgFormat)
%helperConstellationSymbols Get the constellation used in a transmission
%
%   CONST = helperConstellationSymbols(CFGFORMAT) returns the constellation
%   used in a transmission as a column vector of complex symbols. The
%   returned symbols are Gray-coded ordered.
%
%   CFGFORMAT is the format configuration object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a>, 
%   <a href="matlab:help('wlanHTConfig')">wlanHTConfig</a>, or <a href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a>, which specifies the parameters for
%   the VHT, HT-Mixed, and Non-HT formats, respectively.

% Copyright 2015 The MathWorks, Inc.

validateattributes(cfgFormat,{'wlanVHTConfig','wlanHTConfig', ...
    'wlanNonHTConfig'},{'scalar'},mfilename,'format configuration object');

% Determine the modulation and reference from MCS
if isa(cfgFormat,'wlanNonHTConfig')
    coder.internal.errorIf(~strcmp(cfgFormat.Modulation,'OFDM'), ...
        'wlan:helperConstellationSymbols:InvalidNonHTModulation')
    switch cfgFormat.MCS
        case {0,1}   % 'BPSK';  
            N = 2;
            const = pskmod((0:N-1).',N);
        case {2,3}   % 'QPSK';  
            N = 4;
            const = qammod((0:N-1).',N);
        case {4,5}   % '16QAM';
            N = 16;  
            const = qammod((0:N-1).',N);
        otherwise % MCS = {6,7}, '64QAM';   
            N = 64;  
            const = qammod((0:N-1).',N);
    end
else
    if isa(cfgFormat,'wlanVHTConfig')
        modNum = 10; % For VHT max MCS is 9
    else
        modNum = 8;  % For HT modulation scheme mod(MCS,8)
    end

    switch mod(cfgFormat.MCS,modNum)
        case 0       % 'BPSK';  
            N = 2;
            const = pskmod((0:N-1).',N);
        case {1,2}   % 'QPSK';  
            N = 4;
            const = qammod((0:N-1).',N);
        case {3,4}   % '16QAM';
            N = 16;  
            const = qammod((0:N-1).',N);
        case {5,6,7} % '64QAM';   
            N = 64;  
            const = qammod((0:N-1).',N);
        otherwise %  MCS = {8,9}; '256QAM';  
            N = 256;  
            const = qammod((0:N-1).',N);
    end
end

% Normalize the reference constellation
scale = modnorm(const,'AVPOW',1);
const = const*scale;
end