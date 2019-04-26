classdef wlanHTConfig < wlanConfigBase
%wlanHTConfig Create a high throughput (HT) format configuration object
%   CFGHT = wlanHTConfig creates a high throughput (HT) format
%   configuration object. This object contains the transmit parameters for
%   the HT-Mixed Format of the IEEE 802.11 standard.
%
%   CFGHT = wlanHTConfig(Name,Value) creates a HT object, CFGHT, with the
%   specified property Name set to the specified Value. You can specify
%   additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%
%   wlanHTConfig properties:
%
%   ChannelBandwidth     - Channel bandwidth (MHz) 
%   NumTransmitAntennas  - Number of transmit antennas
%   NumSpaceTimeStreams  - Number of space-time streams 
%   NumExtensionStreams  - Number of extension spatial streams
%   SpatialMapping       - Spatial mapping scheme
%   SpatialMappingMatrix - Spatial mapping matrices
%   MCS                  - Modulation and coding scheme
%   GuardInterval        - Guard interval type used for transmission
%   ChannelCoding        - Forward error correction coding type used
%   PSDULength           - Length of the PSDU in bytes
%   RecommendSmoothing   - Recommend smoothing for channel estimation
%
%   % Example: 
%   %  Create a wlanHTConfig object for BCC, SISO operation for a PSDU
%   %  length of 2048 bytes.
%   
%   cfgHT = wlanHTConfig('ChannelBandwidth', 'CBW20');
%   cfgHT.PSDULength = 2048;
%   disp(cfgHT)
%
%   See also wlanWaveformGenerator, wlanVHTConfig, wlanNonHTConfig.

%   Copyright 2015 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

    %%
    % Define properties in order of intended display
    % Public properties
    properties (SetAccess = 'public')
        %ChannelBandwidth Channel bandwidth (MHz) for PPDU transmission
        %   Specify the channel bandwidth for the packet as one of 'CBW20'
        %   | 'CBW40' to indicate 20MHz and 40MHz use respectively. The
        %   default is 'CBW20'.
        ChannelBandwidth = 'CBW20';     
        %NumTransmitAntennas Number of transmit antennas
        %   Specify the number of transmit antennas as a positive integer
        %   scalar between 1 and 4, inclusive.  The default is 1.
        NumTransmitAntennas = 1;        
        %NumSpaceTimeStreams Number of space-time streams
        %   Specify the number of space-time streams as a positive integer
        %   scalar between 1 and 4, inclusive. The default is 1.
        NumSpaceTimeStreams = 1;
        %NumExtensionStreams Number of extension spatial streams
        %   Specify the number of spatial extension streams as a positive
        %   integer scalar between 0 and 3, inclusive. The default is 0.
        NumExtensionStreams = 0; 
        %SpatialMapping Spatial mapping scheme
        %   Specify the spatial mapping scheme as one of 'Direct' |
        %   'Hadamard' | 'Fourier' | 'Custom'. The default value of this
        %   property is 'Direct', which applies when the
        %   NumSpaceTimeStreams and NumTransmitAntennas properties are
        %   equal and the NumExtensionStreams property is 0.
        SpatialMapping = 'Direct';
        %SpatialMappingMatrix Spatial mapping matrices
        %   Specify the spatial mapping matrix(ces) as a double precision,
        %   real or complex, 2D matrix or 3D array. This property applies
        %   when you set the SpatialMapping property to 'Custom'. It can be
        %   of size [Nsts+Ness, Nt], where Nsts is the NumSpaceTimeStreams
        %   property value, Ness is the NumExtensionStreams property value
        %   and Nt is the NumTransmitAntennas property value. In this case,
        %   the spatial mapping matrix applies to all the frequency
        %   subcarriers and its first Nsts and last Ness rows apply to the
        %   space-time streams and extension spatial streams respectively.
        %   Alternatively, it can be of size [Nst, Nsts+Ness, Nt], where
        %   Nst is the number of data plus pilot subcarriers determined by
        %   the ChannelBandwidth property. Specifically, Nst is 56 for
        %   'CBW20' and 114 for 'CBW40'. In this case, each data and pilot
        %   subcarrier can have its own spatial mapping matrix. In either
        %   2D or 3D case, the spatial mapping matrix for each subcarrier
        %   is normalized. The default value of this property is 1.
        SpatialMappingMatrix = 1;
        %MCS Modulation and coding scheme
        %   Specify the modulation and coding scheme for the packet
        %   transmission as a integer scalar between 0 and 31, inclusive.
        %   The selected value also sets the number of spatial streams
        %   (Nss) for the configuration. The difference between the number
        %   of space-time streams (NumSpaceTimeStreams) and Nss conveys
        %   the use of space-time block coding (STBC). The default is 0.
        MCS = 0;                       
        %GuardInterval Guard interval (cyclic prefix) type
        %   Specify the cyclic prefix type of the data field within a
        %   packet as one of 'Long' | 'Short'. An interval of 800ns and
        %   400ns is used for long and short guard interval types
        %   respectively. The default is 'Long'.
        GuardInterval = 'Long';  
    end

    properties (Constant)
        %ChannelCoding Forward error correction coding type
        %   The channel coding type for the data field is set to 'BCC' to 
        %   indicate binary convolutional coding. LDPC coding is not
        %   supported as yet. As a result, this is a read-only property.
        ChannelCoding = 'BCC';      
    end

    % Public properties
    properties (SetAccess = 'public')
        %PSDULength PSDU length
        %   Specify the PSDU length as an integer scalar for the data
        %   carried in a packet, in bytes. The default is 1024.
        PSDULength = 1024;             
        %RecommendSmoothing Recommend smoothing for channel estimation
        %   Set this property to true to indicate smoothing is recommended
        %   for channel estimation. The default is true.
        RecommendSmoothing = true;      
    end

    properties(Constant, Hidden)
        ChannelBandwidth_Values  = {'CBW20', 'CBW40'};
        SpatialMapping_Values    = {'Direct', 'Hadamard', 'Fourier', 'Custom'}
        ChannelCoding_Values     = {'BCC', 'LDPC'};
        GuardInterval_Values     = {'Long', 'Short'};
    end
    
    %%
    methods
        % Constructor
        function obj = wlanHTConfig(varargin)
            % Add sets for enum properties with different sized values
            obj = obj@wlanConfigBase('SpatialMapping','Direct', ...
                                     'GuardInterval', 'Long', varargin{:});
        end
        
        % Property self-validation and sets
        function obj = set.ChannelBandwidth(obj,val)
            prop = 'ChannelBandwidth';
            validateEnumProperties(obj, prop, val);
            obj.(prop) = val;    
        end
               
        function obj = set.NumTransmitAntennas(obj, val)
            prop = 'NumTransmitAntennas';
            validateattributes(val, {'numeric'}, ...
                {'real','integer','scalar' '>=',1,'<=',4}, ...
                [class(obj) '.' prop], prop);
            obj.(prop) = val;
        end

        function obj = set.NumSpaceTimeStreams(obj, val)
            prop = 'NumSpaceTimeStreams';
            validateattributes(val, {'numeric'}, ...
                {'real','integer','scalar' '>=',1,'<=',4}, ...
                [class(obj) '.' prop], prop);
            obj.(prop) = val;
        end

        function obj = set.NumExtensionStreams(obj, val)
            prop = 'NumExtensionStreams';
            validateattributes(val, {'numeric'}, ...
                {'real','integer','scalar','>=',0,'<=',3}, ...
                [class(obj) '.' prop], prop);
            obj.(prop) = val;
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
            coder.internal.errorIf( ...
                (is3DFormat && ~any(size(val, 1) == [56 114])) || ...  
                (numSTS > 4) || (numTx > 4) || (numSTS > numTx), ...
                'wlan:wlanHTConfig:InvalidSpatialMapMtxDim');
            
            obj.(propName) = val;
        end
        
        function obj = set.MCS(obj,val)
            prop = 'MCS';
            validateattributes(val, {'numeric'}, ...
                {'real','integer','scalar','>=',0,'<=',31}, ...
                [class(obj) '.' prop], prop); 
            obj.(prop) = val;
        end
                
        function obj = set.GuardInterval(obj,val)
            prop = 'GuardInterval';
            validateEnumProperties(obj, prop, val);
            obj.(prop) = '';
            obj.(prop) = val;
        end
                                         
        function obj = set.PSDULength(obj,val)
            prop = 'PSDULength';
            validateattributes(val, {'numeric'}, ...
                {'real','integer','scalar','>=',0,'<=',65535}, ...
                [class(obj) '.' prop], prop); 
            obj.(prop) = val;
        end

        function obj = set.RecommendSmoothing(obj,val)
            prop = 'RecommendSmoothing';
            validateattributes(val, {'logical'}, {'scalar'}, ...
                [class(obj) '.' prop], prop);
            obj.(prop) = val;
        end
        
        function varargout = validateConfig(obj, varargin)
            % validateConfig Validate the wlanHTConfig object
            %
            %   validateConfig(CFGHT) validates the dependent properties
            %   for the specified wlanHTConfig configuration object.
            
            %   For INTERNAL use only, subject to future changes:
            %
            %   validateConfig(CFGHT, MODE) validates only the subset of 
            %   dependent properties as specified by the MODE input.
            %   MODE must be one of:
            %    'EssSTS':
            %    'STSTx':
            %    'SMapping':
            %    'MCS':
            %    'MCSSTSTx':
            %    'SMappingMCS':
            
            narginchk(1,2);nargoutchk(0,1);
            if (nargin==2)
                mode = varargin{1};
            else
                mode = 'Full';
            end
            
            s = struct('NumDataSymbols', nan, ...
                'NumPadBits',     nan, ...
                'NumPPDUSamples', nan, ...
                'TxTime',         nan);

            if strcmp(mode, 'EssSTS')           % wlanHTLTFDemodulate, 
                validateEssSTS(obj);            % wlanHTLTFChannelEstimate 
                
            elseif strcmp(mode, 'STSTx')        % wlanFieldIndices
                validateSTSTx(obj);

            elseif strcmp(mode, 'SMapping')     % HT-STF, HT-LTF
                validateSpatialMapping(obj);
                
            elseif strcmp(mode, 'MCS')          % wlanHTDataRecover
                s = validateMCSLength(obj);
                
            elseif strcmp(mode, 'MCSSTSTx')     % HT-SIG, LSIG
                validateSTSTx(obj);
                s = validateMCSLength(obj);
                
            elseif strcmp(mode, 'SMappingMCS') || strcmp(mode, 'Full') % HT-Data
                % Shared full object validation
                validateSpatialMapping(obj);
                s = validateMCSLength(obj);
                               
            end
            if nargout==1
                varargout{1} = s;
            end
                        
        end      
  
    end
    
    methods (Access = protected)
        function flag = isInactiveProperty(obj, prop)
            flag = false;
            if strcmp(prop, 'NumExtensionStreams')
                flag = ~inESSMode(obj);
            elseif strcmp(prop, 'SpatialMappingMatrix')
                flag = ~strcmp(obj.SpatialMapping, 'Custom');
            end
        end
    end

    methods (Access = private)

            function s = privInfo(obj)
            %privInfo Returns information relevant to the object
            %   S = privInfo(cfgHT) returns a structure, S, containing the
            %   relevant information for the wlanHTConfig object, cfgHT.
            %   The output structure S has the following fields:
            %
            %   NumDataSymbols - Number of OFDM symbols for the Data field
            %   NumPadBits     - Number of pad bits in the Data field
            %   NumPPDUSamples - Number of PPDU samples per transmit antenna
            %   TxTime         - PPDU transmission time in us
                        
            % Compute the number of OFDM symbols in Data field
            mcsTable  = getRateTable(obj);
            numDBPS = mcsTable.NDBPS;
            numES   = mcsTable.NES;
            numSS   = mcsTable.Nss;
            
            numSTS = obj.NumSpaceTimeStreams;
            STBC = numSTS - numSS;
            mSTBC = 1 + (STBC~=0);
            
            % Only for Binary convolutional coding for now
            if obj.PSDULength > 0                
                Ntail = 6;
                numDataSym = mSTBC * ceil((8*obj.PSDULength + 16 + ...
                                          Ntail*numES)/(mSTBC*numDBPS));
                numPadBits = numDataSym * numDBPS - (8*obj.PSDULength + ...
                                                     16 + Ntail*numES);
            else % == 0, NDP or sounding packet
                numDataSym = 0;
                numPadBits = 0;
            end

            % Compute the number of PPDU samples at CBW
            switch obj.ChannelBandwidth
                case 'CBW40'
                    Nfft = 128;
                otherwise   % 'CBW20'
                    Nfft = 64;
            end

            NHTDLTFVec = [1 2 4 4]; % only HTDLTFs
            NHTELTFVec = [0 1 2 4]; % only HTELTFS
            if inESSMode(obj)
                numESS = obj.NumExtensionStreams;            
            else
                numESS = 0;
            end
            numPreambSym = 2 + 2 + 1 + 2 + 1 + NHTDLTFVec(numSTS) + ...
                           NHTELTFVec(numESS+1); 
            numSymbols = numPreambSym + numDataSym;
            
            if strcmp(obj.GuardInterval, 'Short')
                cpLen = Nfft/8;
                numSamples = numPreambSym*(Nfft*5/4) + numDataSym*(Nfft + cpLen);
                txTime = numPreambSym*4 + 4*ceil(numDataSym*3.6/4);
            else % 'Long'
                cpLen = Nfft/4;
                numSamples = numSymbols * (Nfft + cpLen);
                txTime = (numPreambSym+numDataSym)*4;
            end
            
            % Set output structure
            s.NumDataSymbols = numDataSym;
            s.NumPadBits     = numPadBits;
            s.NumPPDUSamples = numSamples;
            s.TxTime         = txTime;
        end

        function validateEssSTS(obj)
        %   ValidateESSSTS Validate NumExtensionStreams, NumSpaceTimeStreams
        %   properties for wlanHTConfig object
        
            if inESSMode(obj)
                % Nsts + Ness <= 4
                coder.internal.errorIf( obj.NumExtensionStreams + ...
                    obj.NumSpaceTimeStreams > 4, ...
                    'wlan:wlanHTConfig:InvalidNumEss');
            end
        end

        function validateSTSTx(obj)
        %   ValidateSTSTx Validate NumTransmitAntennas, NumSpaceTimeStreams
        %   NumExtensionStreams properties for wlanHTConfig object
        
            if inESSMode(obj)
                % NumTx and Nsts: numTx cannot be less than (Nsts+Ness)
                coder.internal.errorIf( obj.NumTransmitAntennas < ...
                    obj.NumSpaceTimeStreams + obj.NumExtensionStreams, ...
                    'wlan:wlanHTConfig:InvalidNumTxandSTSESS');
            else
                % NumTx and Nsts: numTx cannot be less than Nsts
                coder.internal.errorIf( obj.NumTransmitAntennas < ...
                    obj.NumSpaceTimeStreams, ...
                    'wlan:wlanHTConfig:InvalidNumTxandSTS');
            end
        end

        function validateSpatialMapping(obj)
        %   ValidateSpatialMapping Validate spatial mapping properties for
        %   wlanHTConfig object that include:
        %     NumTransmitAntennas, NumSpaceTimeStreams, SpatialMapping, 
        %     SpatialMappingMatrix, NumExtensionStreams, ChannelBandwidth

            validateSTSTx(obj);        
            validateEssSTS(obj);
            
            %% Validate spatial mapping between STS+ESS and Tx
            if inESSMode(obj)               
                numESS = obj.NumExtensionStreams;
            else
                numESS = 0;
            end

            coder.internal.errorIf(strcmp(obj.SpatialMapping, 'Direct') && ...
            ( (obj.NumSpaceTimeStreams ~= obj.NumTransmitAntennas) ), ...
            'wlan:wlanHTConfig:NumSTSNotEqualNumTxDirectMap');

            % Restrict to Custom for v1.
            coder.internal.errorIf(numESS ~=0 && ...
                ~strcmp(obj.SpatialMapping, 'Custom'), ...
                'wlan:wlanHTConfig:NumESSNotCustomMap');

            if strcmp(obj.SpatialMapping, 'Custom')
                % Validate spatial mapping matrix
                SPMtx = obj.SpatialMappingMatrix;
                is3DFormat = (ndims(SPMtx) == 3) || (iscolumn(SPMtx) && ...
                             ~isscalar(SPMtx));
                numSTSPlusESS = size(SPMtx, 1+is3DFormat);
                numTx         = size(SPMtx, 2+is3DFormat);
                if strcmp(obj.ChannelBandwidth, 'CBW20')
                    numST = 56;
                else
                    numST = 114;
                end
                coder.internal.errorIf( ...
                    (is3DFormat && (size(SPMtx, 1) ~= numST)) || ...
                    (numSTSPlusESS ~= (obj.NumSpaceTimeStreams + numESS)) || ...
                    (numTx ~= obj.NumTransmitAntennas), ...
                    'wlan:wlanHTConfig:MappingMtxNotMatchOtherProp', ...
                    obj.NumSpaceTimeStreams, numESS, obj.NumTransmitAntennas, ...
                    numST);
            end
        end

        function s = validateMCSLength(obj)
        %   ValidateMCSLength Validate PSDULength, MCS and Spatial properties
        %   for the wlanHTConfig object that include:
        %     PSDULength, MCS, NumSpaceTimeStreams, ChannelBandwidth, 
        %     NumExtensionStreamsm, GuardInterval
        %   Invokes the info method and returns the output, if valid.
        
            validateEssSTS(obj);

            % IEEE Std 802.11-2012, Table 20-12, Nsts cannot be less than Nss
            Nss = floor(obj.MCS/8)+1;
            coder.internal.errorIf( obj.NumSpaceTimeStreams < Nss, ...
                'wlan:wlanHTConfig:InvalidNumSTSandSSLT', ...
                obj.NumSpaceTimeStreams, Nss);

            % IEEE Std 802.11-2012, Tables 20-12, 20-18, Nsts cannot be > 2*Nss
            coder.internal.errorIf( obj.NumSpaceTimeStreams > 2*Nss, ...
                'wlan:wlanHTConfig:InvalidNumSTSandSSGT', ... 
                obj.NumSpaceTimeStreams, Nss);
            
            % Validate PSDULength for txTime (max 5.484ms for Mixed mode)
            s = privInfo(obj);
            coder.internal.errorIf( s.TxTime > 5484, ...
                'wlan:wlanHTConfig:InvalidPPDUDuration'); 
                        
        end
        
    end
end

% [EOF]