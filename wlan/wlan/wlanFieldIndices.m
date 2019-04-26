function indices = wlanFieldIndices(cfgFormat, varargin)
%WLANFIELDINDICES Generate field indices for WLAN packet
%
%   INDICES = wlanFieldIndices(CFGFORMAT) returns the start and end
%   time-domain sample indices for all fields in a packet relative to the
%   first sample in a packet.
%
%   INDICES is a structure array of field names for the specified
%   configuration and contains the start and end indices of all fields in a
%   packet.
%
%   CFGFORMAT is a format configuration object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a>, 
%   <a href="matlab:help('wlanHTConfig')">wlanHTConfig</a>, or <a href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a>.
%
%   INDICES = wlanFieldIndices(CFGFORMAT, FIELDNAME) returns the start and
%   end time-domain sample indices for the specified FIELDNAME in a packet.
%   INDICES is a row vector of size two representing the start and end
%   sample indices of the specified field.
%
%   FIELDNAME is a string specifying the field of interest and depends on
%   the type of CFGFORMAT. FIELDNAME must be one of 'VHT-SIG-A', 'VHT-STF',
%   'VHT-LTF', 'VHT-SIG-B' or 'VHT-Data' for a CFGFORMAT specified as a 
%   <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a> object. For <a href="matlab:help('wlanHTConfig')">wlanHTConfig</a>, the field of interest must be
%   'HT-SIG', 'HT-STF', 'HT-LTF or 'HT-Data'. FIELDNAME 'NonHT-Data' is 
%   valid for <a href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a> object, only OFDM modulation is supported. 
%   FIELDNAME 'L-STF', 'L-LTF' and 'L-SIG' are common for all formats.
% 
%   Example 1:
%   %   Extract the necessary fields to decode VHT data for a received
%   %   802.11ac VHT transmission.
%
%      cfgVHT = wlanVHTConfig;     % Create packet configuration
%      cfgVHT.MCS = 1;             % Modulation: QPSK, Rate: 1/2
%      cfgVHT.APEPLength = 1024;   % A-MPDU length in bytes
%
%      % Generate transmit waveform
%      txWaveform = wlanWaveformGenerator([1;0;0;1],cfgVHT);
%
%      rxWaveform = awgn(txWaveform, 1, 1);              % Add noise
%
%      % Recover transmitted bits, assuming synchronization
%      ind = wlanFieldIndices(cfgVHT);
%      symVHTLTF  = wlanVHTLTFDemodulate(...
%                       rxWaveform(ind.VHTLTF(1):ind.VHTLTF(2),:),cfgVHT);
%      chanEst    = wlanVHTLTFChannelEstimate(symVHTLTF,cfgVHT);
%      rxBits     = wlanVHTDataRecover(rxWaveform(ind.VHTData(1): ...
%                                      ind.VHTData(2),:),chanEst,1,cfgVHT);
%
%   %   Example 2:
%   %   Extract and demodulate L-LTF field from a NonHT OFDM packet.
%      
%      % Create packet configuration
%      cfgNonHT = wlanNonHTConfig('Modulation','OFDM');
%
%      % Generate transmit waveform
%      txWaveform = wlanWaveformGenerator([1;0;0;1],cfgNonHT);
%      rxWaveform = awgn(txWaveform,1,1);              % Add noise
%      indLLTF  = wlanFieldIndices(cfgNonHT,'L-LTF');
%      y = wlanLLTFDemodulate(rxWaveform(indLLTF(1):indLLTF(2)),cfgNonHT);
%
%   See also wlanVHTDataRecover, wlanHTDataRecover, wlanNonHTDataRecover.
 
%   Copyright 2015 The MathWorks, Inc.

%#codegen

% Validate input is a class object
validateattributes(cfgFormat, {'wlanVHTConfig', ...
    'wlanHTConfig','wlanNonHTConfig'}, ...
    {'scalar'}, mfilename, 'format configuration object');

% DSSS modulation in NonHT format is not supported
coder.internal.errorIf(isa(cfgFormat, 'wlanNonHTConfig') &&...
    strcmp(cfgFormat.Modulation, 'DSSS'), ...
    'wlan:wlanFieldIndices:InvalidModulation');

narginchk(1,2)
if nargin == 1
    
    if isa(cfgFormat, 'wlanVHTConfig')
        [numVHTLTF, symLength, N20] = getVHTParams(cfgFormat);
        indLSTF    = getIndices(cfgFormat,'L-STF',N20);
        indLLTF    = getIndices(cfgFormat,'L-LTF',N20);
        indLSIG    = getIndices(cfgFormat,'L-SIG',N20);
        indVHTSIGA = getIndices(cfgFormat,'VHT-SIG-A',N20);
        indVHTSTF  = getIndices(cfgFormat,'VHT-STF',N20);
        indVHTLTF  = getIndices(cfgFormat,'VHT-LTF',N20,numVHTLTF);
        indVHTSIGB = getIndices(cfgFormat,'VHT-SIG-B',N20,numVHTLTF);
        if  isscalar(cfgFormat.APEPLength) && (cfgFormat.APEPLength == 0)
            indVHTData = zeros(0,2,'uint32');
        else
            indVHTData = getIndices(cfgFormat,'VHT-Data',N20, ...
                numVHTLTF,symLength);
        end
        indices = struct('LSTF',    indLSTF, ...
            'LLTF',    indLLTF, ...
            'LSIG',    indLSIG, ...
            'VHTSIGA', indVHTSIGA, ...
            'VHTSTF',  indVHTSTF, ...
            'VHTLTF',  indVHTLTF, ...
            'VHTSIGB', indVHTSIGB, ...
            'VHTData', indVHTData);
        
    elseif isa(cfgFormat,'wlanHTConfig')
        [numHTLTF, symLength, N20] = getHTParams(cfgFormat);
        indLSTF   = getIndices(cfgFormat,'L-STF',N20);
        indLLTF   = getIndices(cfgFormat,'L-LTF',N20);
        indLSIG   = getIndices(cfgFormat,'L-SIG',N20);
        indHTSIG  = getIndices(cfgFormat,'HT-SIG',N20);
        indHTSTF  = getIndices(cfgFormat,'HT-STF',N20);
        indHTLTF  = getIndices(cfgFormat,'HT-LTF',N20,numHTLTF);
        if (cfgFormat.PSDULength ==0)
            indHTData = zeros(0,2,'uint32');
        else
            indHTData = getIndices(cfgFormat,'HT-Data',N20, ...
                numHTLTF,symLength);
        end
        indices = struct('LSTF',  indLSTF, ...
            'LLTF',  indLLTF, ...
            'LSIG',  indLSIG, ...
            'HTSIG', indHTSIG,...
            'HTSTF', indHTSTF,...
            'HTLTF', indHTLTF,...
            'HTData',indHTData);
        
    else
        indLSTF  = getIndices(cfgFormat,'L-STF',1);
        indLLTF  = getIndices(cfgFormat,'L-LTF',1);
        indLSIG  = getIndices(cfgFormat,'L-SIG',1);
        indNonHTData = getIndices(cfgFormat,'NonHT-Data',1);
        
        indices = struct('LSTF',    indLSTF, ...
            'LLTF',     indLLTF, ...
            'LSIG',     indLSIG, ...
            'NonHTData',indNonHTData);
    end
    
else
    fieldType = varargin{1};
    
    if isa(cfgFormat, 'wlanVHTConfig')
        coder.internal.errorIf(~ischar(fieldType) || ...
            ~any(strcmp(fieldType, ...
            {'L-STF','L-LTF','L-SIG','VHT-SIG-A','VHT-STF','VHT-LTF', ...
            'VHT-SIG-B','VHT-Data'})), ...
            'wlan:wlanFieldIndices:InvalidFieldTypeVHT');
        if strcmp(fieldType, 'VHT-Data') && ...
           isscalar(cfgFormat.APEPLength) && (cfgFormat.APEPLength == 0)
            indices = zeros(0,2,'uint32'); % NDP
        else
            [numVHTLTF, symLength, N20] = getVHTParams(cfgFormat);
            indices = getIndices(cfgFormat, fieldType, N20, numVHTLTF, ...
                symLength);
        end
    elseif isa(cfgFormat,'wlanHTConfig')
        coder.internal.errorIf(~ischar(fieldType) || ...
            ~any(strcmp(fieldType, ...
            {'L-STF','L-LTF','L-SIG','HT-SIG','HT-STF','HT-LTF', ...
            'HT-Data'})), ...
            'wlan:wlanFieldIndices:InvalidFieldTypeHT');
        if strcmp(fieldType, 'HT-Data') && (cfgFormat.PSDULength == 0)
            indices = zeros(0,2,'uint32'); % NDP
        else
            [numHTLTF, symLength, N20] = getHTParams(cfgFormat);
            indices = getIndices(cfgFormat, fieldType ,N20, numHTLTF, ...
                symLength);
        end
    else
        
        coder.internal.errorIf(~ischar(fieldType) || ...
            ~any(strcmp(fieldType, ...
            {'L-STF','L-LTF','L-SIG', 'NonHT-Data'})), ...
            'wlan:wlanFieldIndices:InvalidFieldTypeNonHT');
        indices = getIndices(cfgFormat, fieldType, 1);
    end
    
end
end

function out = getIndices(format, fieldType, N20, varargin)

% Start and end sample indices are relative to the field length in samples
% for 20MHz bandwith.

if strcmp(fieldType, 'L-STF')           
    indStart = 1;                 % Start of L-STF field
    indEnd = 160*N20;             % End of L-STF field
elseif strcmp(fieldType, 'L-LTF')     
    indStart = 160*N20+1;         % Start of L-LTF field
    indEnd = 320*N20;             % End of L-LTF field
elseif strcmp(fieldType, 'L-SIG')
    indStart = 320*N20+1;         % Start of L-SIG field
    indEnd = indStart+80*N20-1;   % End of L-SIG field
elseif strcmp(fieldType, 'VHT-SIG-A')
    indStart = 400*N20+1;         % Start of VHT-SIG-A field
    indEnd = indStart+160*N20-1;  % End of VHT-SIG-A field
elseif strcmp(fieldType, 'HT-SIG') 
    indStart = 400*N20+1;         % Start of HT-SIG field
    indEnd = indStart+160*N20-1;  % End of HT-SIG field
elseif strcmp(fieldType, 'VHT-STF')     
    indStart = 560*N20+1;         % Start of VHT-STF field
    indEnd = indStart+80*N20-1;   % End of VHT-STF field
elseif strcmp(fieldType, 'HT-STF')      
    indStart = 560*N20+1;         % Start of HT-STF field
    indEnd = indStart+80*N20-1;   % End of HT-STF field
elseif strcmp(fieldType, 'VHT-LTF')
    numVHTLTF = varargin{1};     
    indStart = 640*N20+1;                 % Start of VHT-LTF field
    indEnd = indStart+numVHTLTF*80*N20-1; % End of VHT-LTF field
elseif strcmp(fieldType, 'HT-LTF')
    numHTLTF = varargin{1}; 
    validateConfig(format, 'EssSTS');
    validateConfig(format, 'STSTx');
    indStart = 640*N20+1;                 % Start of HT-LTF field
    indEnd = indStart+numHTLTF*80*N20-1;  % End of HT-LTF field
elseif strcmp(fieldType, 'VHT-SIG-B')
    numVHTLTF = varargin{1};
    indStart = (640+numVHTLTF*80)*N20+1;   % Start of VHT-SIG-B field
    indEnd = indStart+80*N20-1;           % End of VHT-SIG-B field
elseif strcmp(fieldType, 'VHT-Data')
    S = validateConfig(format, 'MCS'); 
    numVHTLTF = varargin{1};
    symLength = varargin{2};
    indStart = (720+numVHTLTF*80)*N20+1;          % Start of VHT-Data field
    indEnd = indStart+S.NumDataSymbols*symLength-1; % End of VHT-Data field
elseif strcmp(fieldType, 'HT-Data')
    numHTLTF = varargin{1};
    symLength = varargin{2};
    validateConfig(format, 'EssSTS');  
    validateConfig(format, 'MCSSTSTx');
    S = validateConfig(format, 'MCS');
    indStart = (640+numHTLTF*80)*N20+1;           % Start of HT-Data field
    indEnd = indStart+S.NumDataSymbols*symLength-1; % End of HT-Data field
elseif strcmp(fieldType, 'NonHT-Data')
    symLength = 80;
    S = validateConfig(format, 'Full');
    indStart = 400*N20+1; % Start of NonHT-Data field
    % End of NonHT-Data field
    indEnd = indStart+S.NumDataSymbols*symLength-1;
end
out = uint32([indStart, indEnd]);

end

function [numLTF, symLength, N20MHz] = getVHTParams(format)

NVHTLTFVec = [1 2 4 4 6 6 8 8]; % Number of VHT-LTF symbols
numLTF = NVHTLTFVec(sum(format.NumSpaceTimeStreams));
channelBandwidth = format.ChannelBandwidth;

N20MHz = real(str2double(channelBandwidth(4:end)))/20;
FFTLen = 64 * N20MHz; % FFT length for 20MHz bandwidth

if strcmp(format.GuardInterval,'Short')
    % Short symbol length for 20MHz bandwidth
    symLength = FFTLen * 9/8; 
else
    % Long symbol length for 20MHz bandwidth
    symLength = FFTLen * 5/4;
end

end
    
function [numLTF, symLength, N20MHz] = getHTParams(format)

chanBW = format.ChannelBandwidth;
if inESSMode(format)
    numESS = format.NumExtensionStreams;
else
    numESS = 0;
end

[~, ~, Ndltf, Neltf] = vhtltfSequence(chanBW, ...
    format.NumSpaceTimeStreams, numESS);
numLTF = Ndltf + Neltf;

N20MHz = real(str2double(chanBW(4:end)))/20;
FFTLen = 64 * N20MHz; % FFT length for 20MHz bandwidth

if strcmp(format.GuardInterval,'Short')
    % Short symbol length for 20MHz bandwidth
    symLength = FFTLen * 9/8;
else
    % Long symbol length for 20MHz bandwidth
    symLength = FFTLen * 5/4;
end

end

% [EOF]
