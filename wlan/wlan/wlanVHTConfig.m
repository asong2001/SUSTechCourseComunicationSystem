classdef wlanVHTConfig < wlanConfigBase 
%wlanVHTConfig Create a very high throughput (VHT) format configuration object
%   CFGVHT  = wlanVHTConfig creates a very high throughput format
%   configuration object. This object contains the transmit parameters for
%   the VHT format of IEEE 802.11 standard.
%
%   CFGVHT = wlanVHTConfig(Name,Value) creates a VHT object, CFGVHT, with
%   the specified property Name set to the specified value. You can specify
%   additional name-value pair arguments in any order as (Name1,Value1,
%   ...,NameN,ValueN).
%
%   wlanVHTConfig properties:
%
%   ChannelBandwidth     - Channel bandwidth 
%   NumUsers             - Number of users
%   UserPositions        - User positions
%   NumTransmitAntennas  - Number of transmit antennas 
%   NumSpaceTimeStreams  - Number of space-time streams 
%   SpatialMapping       - Spatial mapping scheme
%   SpatialMappingMatrix - Spatial mapping matrix(ces)
%   Beamforming          - Enable beamforming
%   STBC                 - Enable space-time block coding
%   MCS                  - Modulation and coding schemes
%   ChannelCoding        - Channel coding
%   APEPLength           - APEP lengths
%   PSDULength           - Number of bytes to be coded in the packet 
%                          including the A-MPDU and any MAC padding
%   GuardInterval        - Guard interval type
%   GroupID              - Group identifier
%   PartialAID           - Partial association identifier 
%
%   % Example 1: 
%   %  Create wlanVHTConfig object for a 40MHz, single-user configuration
%
%   cfgVHT = wlanVHTConfig;
%   cfgVHT.ChannelBandwidth = 'CBW40';
%   disp(cfgVHT) 
%
%   % Example 2: 
%   %  Create wlanVHTConfig object for a 20MHz, two-user configuration.
%   %  Each element of vector-valued properties apply to a specific user.
%
%   cfgMU = wlanVHTConfig('ChannelBandwidth', 'CBW20', ...
%                         'NumUsers', 2, ...
%                         'GroupID', 2, ...
%                         'NumTransmitAntennas', 2);
%   cfgMU.NumSpaceTimeStreams = [1 1];
%   cfgMU.MCS                 = [4 8];
%   cfgMU.APEPLength          = [1024 2048];
%   disp(cfgMU) 
%
%   See also wlanHTConfig, wlanNonHTConfig.

%   Copyright 2015 The MathWorks, Inc.

%#codegen
    
properties (Access = 'public')
    %ChannelBandwidth Channel bandwidth (MHz) of PPDU transmission
    %   Specify the channel bandwidth as one of 'CBW20' | 'CBW40' | 'CBW80'
    %   | 'CBW160'. The default value of this property is 'CBW80'.
    ChannelBandwidth = 'CBW80';
    %NumUsers Number of users
    %   Specify the number of users as an integer scalar between 1 and 4,
    %   inclusive. The default value of this property is 1.
    NumUsers = 1;
    %UserPositions User positions
    %   Specify the user positions as an integer row vector with length
    %   equal to NumUsers and elements between 0 and 3, inclusive, in a
    %   strictly increasing order. This property applies when you set the
    %   NumUsers property to 2, 3 or 4. The default value of this property
    %   is [0 1].
    UserPositions = [0 1];
    %NumTransmitAntennas Number of transmit antennas
    %   Specify the number of transmit antennas as an integer scalar
    %   between 1 and 8, inclusive. The default value of this property is
    %   1.
    NumTransmitAntennas = 1;
    %NumSpaceTimeStreams Number of space-time streams per user
    %   Specify the number of space-time streams as integer scalar or row
    %   vector with length equal to NumUsers. For a scalar, it must be
    %   between 1 and 8, inclusive. For a row vector, all elements must be
    %   between 1 and 4, inclusive, and sum to no larger than 8. The
    %   default value of this property is 1.
    NumSpaceTimeStreams = 1;
    %SpatialMapping Spatial mapping scheme
    %   Specify the spatial mapping scheme as one of 'Direct' | 'Hadamard'
    %   | 'Fourier' | 'Custom'. The default value of this property is
    %   'Direct', which applies when the sum of the elements in
    %   NumSpaceTimeStreams is equal to NumTransmitAntennas.
    SpatialMapping = 'Direct';
    %SpatialMappingMatrix Spatial mapping matrix(ces)
    %   Specify the spatial mapping matrix(ces) as a real or complex, 2D
    %   matrix or 3D array. This property applies when you set the
    %   SpatialMapping property to 'Custom'. It can be of size
    %   NstsTotal-by-Nt, where NstsTotal is the sum of the elements in the
    %   NumSpaceTimeStreams property and Nt is the NumTransmitAntennas
    %   property. In this case, the spatial mapping matrix applies to all
    %   the subcarriers. Alternatively, it can be of size Nst-by-
    %   NstsTotal-Nt, where Nst is the number of occupied subcarriers
    %   determined by the ChannelBandwidth property. Specifically, Nst is
    %   56 for 'CBW20', 114 for 'CBW40', 242 for 'CBW80' and 484 for
    %   'CBW160'. In this case, each occupied subcarrier can have its own
    %   spatial mapping matrix. In either 2D or 3D case, the spatial
    %   mapping matrix for each subcarrier is normalized. The default value
    %   of this property is 1.
    SpatialMappingMatrix = 1;
    %Beamforming Enable beamforming
    %   Set this property to true when the specified SpatialMappingMatrix
    %   property is a beamforming steering matrix(ces). This property
    %   applies when you set the NumUsers property to 1 and the
    %   SpatialMapping property to 'Custom'. The default value of this
    %   property is true.
    Beamforming = true;
    %STBC Enable space-time block coding
    %   Set this property to true to enable space-time block coding in the
    %   data field transmission. This property applies when you set the
    %   NumUsers property to 1. The default value of this property is
    %   false.
    STBC = false;
    %MCS Modulation and coding scheme per user 
    %   Specify the modulation and coding scheme per user as an integer
    %   scalar or row vector with length equal to NumUsers. Its elements
    %   must be integers between 0 and 9, inclusive. A scalar value applies
    %   to all users. The default value of this property is 0.
    MCS = 0;
end

properties (SetAccess = private)
    %ChannelCoding Forward error correction code type
    %   This is a read-only property. Only binary convolutional coding
    %   (BCC) is supported.
    ChannelCoding = 'BCC';
end 

properties (SetAccess = 'public')  
    %APEPLength APEP length per user
    %   Specify the APEP length in bytes per user as an integer scalar or
    %   row vector with length equal to NumUsers. A scalar value applies to
    %   all users. All elements must be integers between 1 and 1048575,
    %   inclusive. In addition, this property can be 0 when the NumUsers
    %   property is 1, which implies a VHT non-data-packet (NDP). The
    %   default value of this property is 1024.
    APEPLength = 1024;
end

properties (SetAccess = private, GetAccess = public)
    %PSDULength PSDU lengths
    %   The number of bytes carried in a packet, including the A-MPDU and
    %   any MAC padding. This property is read-only and is calculated
    %   internally based on other properties.
    PSDULength;
end

properties (Access = 'public') 
    %GuardInterval Guard interval type
    %   Specify the guard interval (cyclic prefix) type for data field
    %   transmission as one of 'Long' | 'Short'. The default value of this
    %   property is 'Long'.
    GuardInterval = 'Long';
    %GroupID Group identifier
    %   Specify the group identifier as an integer scalar. It must be 0 or
    %   63 when the NumUsers property is 1 and must be between 1 and 62
    %   inclusive when the NumUsers property is 2, 3 or 4. The default
    %   value of this property is 63.
    GroupID = 63;
    %PartialAID Partial association identifier 
    %   Specify the partial association identifier of the intended
    %   recipient as an integer scalar between 0 and 511, inclusive. This
    %   property applies when you set the NumUsers property to 1. For an
    %   uplink transmission, it is the last nine bits of the BSSID. For a
    %   downlink transmission, it combines the association ID and the BSSID
    %   of its serving AP. The default value of this property is 275.
    PartialAID = 275;
end

properties(Constant, Hidden)
    ChannelBandwidth_Values  = {'CBW20', 'CBW40','CBW80','CBW160'};
    SpatialMapping_Values    = {'Direct', 'Hadamard', 'Fourier', 'Custom'}
    GuardInterval_Values     = {'Short', 'Long'};
end

methods
  function obj = wlanVHTConfig(varargin)
    obj@wlanConfigBase( ...
        'ChannelBandwidth', 'CBW80',  ...
        'SpatialMapping',   'Direct', ...
        'GuardInterval',    'Long',   ...
        varargin{:});
  end

  function obj = set.ChannelBandwidth(obj,val)
    propName = 'ChannelBandwidth';
    validateEnumProperties(obj, propName, val);
    obj.(propName) = ''; 
    obj.(propName) = val;
  end

  function obj = set.NumUsers(obj, val)
    propName = 'NumUsers';
    validateattributes(val, {'numeric'}, ...
        {'real','integer','scalar','>=',1,'<=',4}, ...
        [class(obj) '.' propName], propName); 
    obj.(propName)= val;
  end

  function obj = set.UserPositions(obj, val)
    propName = 'UserPositions';
    validateattributes(val, {'numeric'}, ...
        {'real','integer','row','>=',0,'<=',3,'increasing'}, ...
        [class(obj) '.' propName], propName);
    obj.(propName) = val;                
  end

  function obj = set.NumTransmitAntennas(obj, val)
    propName = 'NumTransmitAntennas';
    validateattributes(val, {'numeric'}, ...
        {'real','integer','scalar','>=',1,'<=',8}, ...
        [class(obj) '.' propName], propName);
    obj.(propName) = val;
  end

  function obj = set.NumSpaceTimeStreams(obj, val)
    propName = 'NumSpaceTimeStreams';
    validateattributes(val, {'numeric'}, ...
        {'real','integer','row','>=',1,'<=',8}, ...
        [class(obj) '.' propName], propName);

    coder.internal.errorIf(~isscalar(val) && ...
        ((length(val) > 4) || any(val > 4) || sum(val) > 8), ...
        'wlan:wlanVHTConfig:InvalidMUSTS'); 

    obj.(propName) = val;
  end

  function obj = set.SpatialMapping(obj, val)
    propName = 'SpatialMapping';
    validateEnumProperties(obj, propName, val);
    obj.(propName) = ''; 
    obj.(propName) = val;
  end

  function obj = set.SpatialMappingMatrix(obj, val)
    propName = 'SpatialMappingMatrix';
    validateattributes(val, {'double'}, {'3d','finite','nonempty'}, ...
        [class(obj) '.' propName], propName); 

    is3DFormat = (ndims(val) == 3) || (iscolumn(val) && ~isscalar(val));
    numSTS = size(val, 1+is3DFormat);
    numTx  = size(val, 2+is3DFormat);
    numST = [56 114 242 484]; % Total number of occupied subcarriers
    coder.internal.errorIf( ...
        (is3DFormat && ~any(size(val, 1) == numST)) || ...  
        (numSTS > 8) || (numTx > 8) || (numSTS > numTx), ...
        'wlan:wlanVHTConfig:InvalidSpatialMapMtxDim');

    obj.(propName) = val;
  end
  
  function obj = set.Beamforming(obj, val)
    propName = 'Beamforming';
    validateattributes(val, {'logical'}, {'scalar'}, ...
        [class(obj) '.' propName], propName);
    obj.(propName) = val;
  end

  function obj = set.STBC(obj, val)
    propName = 'STBC';
    validateattributes(val, {'logical'}, {'scalar'}, ...
        [class(obj) '.' propName], propName);
    obj.(propName) = val;
  end

  function obj = set.MCS(obj, val)
    propName = 'MCS';
    validateattributes(val, {'numeric'}, ...
      {'real','integer','row','>=',0,'<=',9}, ...
      [class(obj) '.' propName], propName);

    coder.internal.errorIf(length(val) > 4, ...
      'wlan:wlanVHTConfig:InvalidMUMCS');

    obj.(propName) = val;
  end

  function obj = set.APEPLength(obj, val)
    propName = 'APEPLength';
    maxBytes = 1048575; % Maximum number of bytes
    validateattributes(val, {'numeric'}, ...
      {'real','integer','row','>=',0,'<=',maxBytes}, ...
      [class(obj) '.' propName], propName);

    coder.internal.errorIf(~isscalar(val) && ...
      ((length(val) > 4) || any(val == 0)), ...
      'wlan:wlanVHTConfig:InvalidMUAPEPLen');

    obj.(propName) = val;
  end

  function obj = set.GuardInterval(obj, val)
    propName = 'GuardInterval';
    validateEnumProperties(obj, propName, val);
    obj.(propName) = ''; 
    obj.(propName) = val;
  end

  function obj = set.GroupID(obj, val)
    propName = 'GroupID';
    validateattributes(val, {'numeric'}, ...
        {'real','scalar','integer','>=',0,'<=',63}, ...
        [class(obj) '.' propName], propName);
    obj.(propName) = val;
  end

  function obj = set.PartialAID(obj, val)
    propName = 'PartialAID';
    validateattributes(val, {'numeric'}, ...
        {'real','scalar','integer','>=',0,'<=',511}, ...
        [class(obj) '.' propName], propName);
    obj.(propName) = val;
  end

  function PSDULen = get.PSDULength(obj)
    % Returns PSDU length in bytes for all users.
    
    if isscalar(obj.APEPLength) && (obj.APEPLength ==0)
        PSDULen = 0;
    else
        validateConfig(obj,'MCSSTSTx');

        numUsers = obj.NumUsers;
        APEPLen  = repmat(obj.APEPLength, 1, numUsers/length(obj.APEPLength));
        mcsTable = getRateTable(obj);
        numDBPS  = mcsTable.NDBPS;
        numES    = mcsTable.NES;

        % Calculate number of OFDM symbols
        numTailBits = 6;
        mSTBC = (numUsers == 1) * (obj.STBC ~= 0) + 1;
        numDataSymbols = max(mSTBC * ceil((8*APEPLen + 16 + ...
                             numTailBits * numES)./(mSTBC * numDBPS)));
        PSDULen = floor((numDataSymbols * numDBPS - 16 - ...
                         numTailBits * numES)/8);
    end
  end

  function varargout = validateConfig(obj, varargin)
    % validateConfig Validate the wlanVHTConfig object
    %   validateConfig(CFGVHT) validates the dependent properties for the
    %   specified wlanVHTConfig configuration object.
    % 
    %   For INTERNAL use only, subject to future changes:
    %
    %   validateConfig(CFGVHT, MODE) validates only the subset of dependent
    %   properties as specified by the MODE input. MODE must be one of:
    %       'SpatialMCSGID'
    %       'SMapping'
    %       'MCS'
    %       'MCSSTSTx'
    %       'SMappingMCS'

    narginchk(1,2);
    nargoutchk(0,1);
    if (nargin==2)
        mode = varargin{1};
    else
        mode = 'Full';
    end

    switch mode
        case 'SpatialMCSGID'    % VHT-SIG-A
            validateSpatial(obj);
            s = validateMCSLength(obj);
            validateGIDAndUserPos(obj);
            
        case 'SMapping'     % VHT-STF, VHT-LTF
            validateSpatialMapping(obj);
            
        case 'MCS'          % wlanVHTDataRecover
            s = validateMCSLength(obj);
            
        case 'MCSSTSTx'     % L-SIG
            validateSTSTx(obj);
            s = validateMCSLength(obj);
            
        case 'SMappingMCS'  % VHT-SIG-B, VHT-Data
            validateSpatialMapping(obj);
            s = validateMCSLength(obj);
            
        otherwise         % wlanWaveformGenerator
            % Full object validation
            validateSpatialMapping(obj);
            s = validateMCSLength(obj);
            validateGIDAndUserPos(obj);
    end

    if nargout == 1
        varargout{1} = s;
    end                        
  end

end   

methods (Access = protected)
  function flag = isInactiveProperty(obj, prop)
    flag = false;
    if any(strcmp(prop, {'STBC', 'PartialAID'}))
        flag = (obj.NumUsers > 1);
    elseif strcmp(prop, 'UserPositions')
        flag = (obj.NumUsers == 1);
    elseif strcmp(prop, 'SpatialMappingMatrix')
        flag = ~strcmp(obj.SpatialMapping, 'Custom');
    elseif strcmp(prop, 'Beamforming')
        flag = (obj.NumUsers > 1) || ~strcmp(obj.SpatialMapping, 'Custom');
    elseif strcmp(prop, 'PSDULength')  
        % Always hide the property
        flag = true;
    end
  end
end

methods (Access = private)
  function s = privInfo(obj)
    %privInfo Returns information relevant to the object
    %   S = privInfo(cfgVHT) returns a structure, S, containing the
    %   relevant information for the wlanVHTConfig object, cfgVHT.
    %   The output structure S has the following fields:
    %
    %   NumDataSymbols - Number of OFDM symbols for the Data field
    %   NumPadBits     - Number of pad bits in the Data field
    %   NumPPDUSamples - Number of PPDU samples per transmit antennas
    %   TxTime         - The time in microseconds, required to
    %                    transmit the PPDU.

    % Set output structure
    numUsers = obj.NumUsers;
    APEPLen  = repmat(obj.APEPLength, 1, numUsers/length(obj.APEPLength));            
    mcsTable = getRateTable(obj);
    numDBPS  = mcsTable.NDBPS;
    numES    = mcsTable.NES;

    % Calculate number of OFDM symbols
    if isscalar(obj.APEPLength) && (obj.APEPLength == 0) % NDP 
        numDataSymbols = 0;
        numPadBits = 0;
    else
        numTailBits = 6;
        mSTBC = (numUsers == 1)*(obj.STBC ~= 0) + 1;
        numDataSymbols = max(mSTBC*ceil((8*APEPLen + 16 + ...
                             numTailBits*numES)./(mSTBC*numDBPS)));
        PSDULen = floor((numDataSymbols*numDBPS - 16 - numTailBits*numES)/8);
        numPadBits = numDataSymbols*numDBPS -  ...
                        (8*PSDULen + 16 + numTailBits*numES);
    end
    
    % Calculate burst time
    NVHTLTFVec = [1 2 4 4 6 6 8 8];
    numPreambSym = 4 + 1 + 2 + 1 + NVHTLTFVec(sum(obj.NumSpaceTimeStreams)) + 1;            
    FFTLen = 64 * str2double(obj.ChannelBandwidth(4:end))/20;
    if strcmp(obj.GuardInterval, 'Short')
        txTime = 4*numPreambSym + 4*ceil(numDataSymbols*3.6/4);
        numPPDUSamples = numPreambSym*FFTLen*5/4 + numDataSymbols*FFTLen*9/8;
    else
        txTime = 4*numPreambSym + 4*numDataSymbols;
        numPPDUSamples = (numPreambSym + numDataSymbols)*FFTLen*5/4;
    end            

    s = struct(...
        'NumDataSymbols', numDataSymbols, ...
        'NumPadBits',     numPadBits, ...
        'NumPPDUSamples', numPPDUSamples, ...
        'TxTime',         txTime); 
  end
        
  function validateSTSTx(obj)
    %   ValidateSTSTx Validate NumTransmitAntennas, NumSpaceTimeStreams
    %   properties for wlanVHTConfig object

    % NumTx and Nsts: numTx cannot be less than sum(Nsts)
    coder.internal.errorIf(obj.NumTransmitAntennas < sum(obj.NumSpaceTimeStreams), ...
        'wlan:wlanVHTConfig:NumSTSLargerThanNumTx');
  end

  function validateSpatial(obj)
    %   validateSpatial Validate the spatial properties for the 
    %   wlanVHTConfig object
    %   Validated property-subset includes:
    %     NumTransmitAntennas, NumSpaceTimeStreams, SpatialMapping

    validateSTSTx(obj);

    coder.internal.errorIf(strcmp(obj.SpatialMapping, 'Direct') && ...
        (sum(obj.NumSpaceTimeStreams) ~= obj.NumTransmitAntennas), ...
        'wlan:wlanVHTConfig:NumSTSNotEqualNumTxDirectMap');            
  end        

  function validateSpatialMapping(obj)
    %   validateSpatialMapping Validate the spatial mapping properties for
    %   the wlanVHTConfig object    
    %   Validated property-subset includes:
    %     ChannelBandwidth, NumTransmitAntennas, NumSpaceTimeStreams, ..
    %     SpatialMapping, SpatialMappingMatrix

    validateSpatial(obj);

    if strcmp(obj.SpatialMapping, 'Custom')
        % Validate spatial mapping matrix
        SPMtx = obj.SpatialMappingMatrix;
        is3DFormat  = (ndims(SPMtx) == 3) || (iscolumn(SPMtx) && ~isscalar(SPMtx));
        numSTSTotal = size(SPMtx, 1+is3DFormat);
        numTx       = size(SPMtx, 2+is3DFormat);
        switch obj.ChannelBandwidth
          case 'CBW20'
            numST = 56;
          case 'CBW40'
            numST = 114;
          case 'CBW80'
            numST = 242;
          otherwise
            numST = 484;
        end
        coder.internal.errorIf( ...
            (is3DFormat && (size(SPMtx, 1) ~= numST)) || ...  
            (numSTSTotal ~= sum(obj.NumSpaceTimeStreams))  || ...
            (numTx ~= obj.NumTransmitAntennas), ...
            'wlan:wlanVHTConfig:MappingMtxNotMatchOtherProp', ...
            sum(obj.NumSpaceTimeStreams), obj.NumTransmitAntennas, numST);
    end            
  end        

  function s = validateMCSLength(obj)
    %ValidateMCSLength Validate MCS and Length properties for
    %   wlanVHTConfig configuration object
    %   Validated property-subset includes:   
    %     ChannelBandwidth, NumUsers, NumSpaceTimeStreams, STBC, MCS,
    %     ChannelCoding, GuardInterval, APEPLength

    numUsers = obj.NumUsers;

    coder.internal.errorIf(length(obj.NumSpaceTimeStreams) ~= ...
        obj.NumUsers, 'wlan:wlanVHTConfig:InvalidSTSNumUsers');

    coder.internal.errorIf((numUsers == 1) && obj.STBC && ...
        all(mod(obj.NumSpaceTimeStreams, 2) == 1), ...
        'wlan:wlanVHTConfig:OddNumSTSWithSTBC');

    coder.internal.errorIf(all(length(obj.MCS) ~= [1 numUsers]), ...
        'wlan:wlanVHTConfig:InvalidMCSNumUsers');

    coder.internal.errorIf(all(length(obj.APEPLength) ~= [1 numUsers]) || ...
        ((numUsers > 1) && any(obj.APEPLength == 0)), ...
        'wlan:wlanVHTConfig:InvalidAPEPLenNumUsers');

    % Check Bandwidth/MCS/Nss valid combinations
    %   Reference: Tables 22-30:22-56, IEEE Std 802.11ac-2013
    invalidComb = ...
          [20,  9, 1; ... % [chanBW, MCS, numSS]
           20,  9, 2; ...
           20,  9, 4; ...
           20,  9, 5; ...
           20,  9, 7; ...
           20,  9, 8; ...
           80,  6, 3; ...
           80,  6, 7; ...
           80,  9, 6; ...
           160, 9, 3];           
    chanBW = str2double(obj.ChannelBandwidth(4:end)); 
    vecMCS = repmat(obj.MCS, 1, numUsers/length(obj.MCS));
    numSS  = obj.NumSpaceTimeStreams / (((numUsers == 1) && obj.STBC) + 1);
    for u = 1:numUsers
        thisComb = [chanBW, vecMCS(u), numSS(u)];
        coder.internal.errorIf(any(ismember(invalidComb, thisComb, 'rows')), ...
            'wlan:wlanVHTConfig:InvalidMCSCombination', ...
            ['''', obj.ChannelBandwidth, ''''], ...
            obj.NumSpaceTimeStreams(u), vecMCS(u), u);
    end

    % Validate PSDULength for txTime (max 5.484ms for VHT format)
    s = privInfo(obj);
    coder.internal.errorIf(s.TxTime > 5484, ...
        'wlan:wlanVHTConfig:InvalidPPDUDuration'); 
  end

  function validateGIDAndUserPos(obj)
    %   validateGIDAndUserPos Validate UserPositions and GroupID against
    %   NumUsers for wlanVHTConfig object.    
    %   Validated property-subset includes:
    %       NumUsers, UserPositions, GroupID

    coder.internal.errorIf((obj.NumUsers > 1) && ...
        (length(obj.UserPositions) ~= obj.NumUsers), ...
        'wlan:wlanVHTConfig:InvalidUserPosNumUsers');

    coder.internal.errorIf(xor((obj.NumUsers == 1), ...
        any(obj.GroupID == [0 63])), ...
        'wlan:wlanVHTConfig:InvalidGIDNumUsers');
  end  
end
    
end

% [EOF]