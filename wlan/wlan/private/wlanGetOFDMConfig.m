function [ofdmCfg, varargout] = wlanGetOFDMConfig(chanBW, CPType, fieldType, varargin)
% WLANGETOFDMCONFIG obtains OFDM configuration parameters.
% 
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   OFDMCFG = wlanGetOFDMConfig(CHANBW, CPTYPE, FIELDTYPE, NUMSTS) returns
%   the OFDM configuration as a structure.
%
%   CHANBW is a string describing the channel bandwidth which must be one
%   of the following: 'CBW20','CBW40','CBW80','CBW160'. 'CBW10', 'CBW5' are
%   also supported and result in same configuration as 'CBW20'.
%
%   CPTYPE is a string and must be one of 'Long','Short'.
%
%   FIELDTYPE is a string and must be one of 'Legacy','HT','VHT'.
%
%   NUMSTS is the number of space-time streams.
%
%   [...,DATAINDNST,PILOTINDNST] = wlanGetOFDMConfig(...) additionally
%   returns the indices of data and pilots within the occupied subcarriers.
%   Both data and pilot indices are column vectors.

% Copyright 2015 The MathWorks, Inc.

%#codegen

narginchk(3,4);
nargoutchk(0,3);

if ~isempty(varargin)
    numSTS = varargin{1}; % Need numSTS for HT or VHT (non-legacy)
else
    numSTS = 1;
end

switch chanBW
  case 'CBW40'
    Num20MHzChan = 2; 
    GammaPhase = [1 1i];                
  case {'CBW80', 'CBW80+80'}
    Num20MHzChan = 4; 
    GammaPhase = [1 -1 -1 -1];
  case 'CBW160'
    Num20MHzChan = 8; 
    GammaPhase = [1 -1 -1 -1 1 -1 -1 -1];
  otherwise   % For 'CBW20', 'CBW10', 'CBW5'
    Num20MHzChan = 1; 
    GammaPhase = 1;
end

FFTLen = 64*Num20MHzChan;
% Carrier rotation: IEEE Std 802.11ac-2013 Section 22.3.7.5.
carrierRotation = reshape(repmat(GammaPhase, [64, 1]), [], 1);

if strcmp(CPType, 'Long')
    CPLen = FFTLen/4;
elseif strcmp(CPType, 'Short')
    CPLen = FFTLen/8;
else % 'Half'
    CPLen = FFTLen/2;
end

if strcmp(fieldType, 'Legacy') 
    numGuardBands = [6; 5]; 
    pilotIdx20MHz =  [12; 26; 40; 54];
    normFactor = FFTLen/sqrt(52*Num20MHzChan);
    
    % Get non-data subcarrier indices per 20MHz channel bandwidth
    nonDataIdxPerGroup = [(1:numGuardBands(1))'; 33; ...
        (64-numGuardBands(2)+1:64)'; pilotIdx20MHz];
    % Get non-data subcarrier indices for the whole bandwidth
    nonDataIdxAll = bsxfun(@plus, nonDataIdxPerGroup, 64*(0:Num20MHzChan-1));
    dataIdx  = setdiff((1:FFTLen)', sort(nonDataIdxAll(:)));
    pilotIdx = reshape(bsxfun(@plus, pilotIdx20MHz, 64*(0:Num20MHzChan-1)), [], 1);
else % 'HT' & 'VHT'
    switch chanBW
      case 'CBW40'
        numGuardBands = [6; 5]; 
        customNullIdx = [-1; 1];
        pilotIdx = [-53; -25; -11; 11; 25; 53];  
      case 'CBW80'
        numGuardBands = [6; 5];
        customNullIdx = [-1; 1];
        pilotIdx = [-103; -75; -39; -11; 11; 39; 75; 103];
      case 'CBW80+80' % Merge with 80MHz case, if same. Separate for now
        numGuardBands = [6; 5];
        customNullIdx = [-1; 1];
        pilotIdx = [-103; -75; -39; -11; 11; 39; 75; 103]; 
      case 'CBW160'
        numGuardBands = [6; 5];
        customNullIdx = [(-129:-127)'; (-5:-1)'; (1:5)'; (127:129)'];
        pilotIdx = [-231; -203; -167; -139; -117; -89; -53; -25; 25; 53; 89; 117; 139; 167; 203; 231]; 
      otherwise  % CBW20, CBW10, CBW5
        numGuardBands = [4; 3];
        customNullIdx = [];
        pilotIdx = [-21; -7; 7; 21];
    end
    
    pilotIdx = pilotIdx + FFTLen/2 + 1; % Convert to 1-based indexing
    customNullIdx = customNullIdx + FFTLen/2 + 1; % Convert to 1-based indexing
    % Get non-data subcarrier indices for the whole bandwidth
    nonDataIdx = [(1:numGuardBands(1))'; FFTLen/2+1; ...
        (FFTLen-numGuardBands(2)+1:FFTLen)'; pilotIdx; customNullIdx];
    dataIdx = setdiff((1:FFTLen)', sort(nonDataIdx));
    normFactor = FFTLen/sqrt(numSTS*(length(dataIdx)+length(pilotIdx)));
end

ofdmCfg = struct( ...
    'FFTLength',           FFTLen, ...
    'CyclicPrefixLength',  CPLen, ...
    'DataIndices',         dataIdx, ...
    'PilotIndices',        pilotIdx, ...
    'CarrierRotations',    carrierRotation, ...
    'NormalizationFactor', normFactor);

if nargout>1
    % Transform indices addressing whole FFT length, to indices addressing
    % occupied subcarriers
    allIndices = [dataIdx; pilotIdx];
    Nsd = numel(dataIdx);
    [~,idxOccupiedSubcarriers] = ismember(allIndices,sort(allIndices));
    dataIndNst = idxOccupiedSubcarriers(1:Nsd); % Data indices within occupied subcarriers
    varargout{1} = dataIndNst;
    if nargout>2
        pilotIndNst = idxOccupiedSubcarriers(Nsd+1:end); % Pilot indices within occupied subcarriers
        varargout{2} = pilotIndNst;
    end
end
end

% [EOF]