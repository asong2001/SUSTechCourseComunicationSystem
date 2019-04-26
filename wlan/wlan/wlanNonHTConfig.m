classdef wlanNonHTConfig < wlanConfigBase
    %wlanNonHTConfig Create a Non-HT format configuration object
    %   CFGNONHT = wlanNonHTConfig creates a Non-HT format configuration
    %   object. This object contains the transmit parameters for the OFDM
    %   and DSSS modulations of the IEEE 802.11 a/b/g standard.
    %
    %   CFGNONHT = wlanNonHTConfig(Name,Value) creates a Non-HT object,
    %   CFGNONHT, with the specified property Name set to the specified
    %   Value. You can specify additional name-value pair arguments in any
    %   order as (Name1,Value1,...,NameN,ValueN).
    %
    %   wlanNonHTConfig properties:
    %
    %   Modulation          - Modulation type
    %   ChannelBandwidth    - Channel bandwidth (MHz) for OFDM modulation
    %   MCS                 - Modulation and coding scheme for OFDM modulation
    %   DataRate            - Date rate for DSSS modulation
    %   Preamble            - Preamble type for DSSS modulation
    %   LockedClocks        - Clock locking indication for DSSS modulation
    %   PSDULength          - Length of the PSDU in bytes
    %   NumTransmitAntennas - Number of transmit antennas
    %
    %   Example: 
    %   %  Create a wlanNonHTConfig object for OFDM operation for a PSDU
    %   %  length of 2048 bytes.
    %
    %   cfgNonHT = wlanNonHTConfig('Modulation', 'OFDM');
    %   cfgNonHT.PSDULength = 2048;
    %   disp(cfgNonHT)
    %
    %   See also wlanWaveformGenerator, wlanVHTConfig, wlanHTConfig.
    
    %   Copyright 2015 The MathWorks, Inc.
        
    %#codegen
    
    % Define properties in order of intended display
    % Public properties
    properties (SetAccess = 'public')
        %Modulation Modulation type 
        %   Specify the modulation type as one of 'OFDM' | 'DSSS' for the
        %   Non-HT transmission packet. The default is 'OFDM'.
        Modulation = 'OFDM';

        %ChannelBandwidth Channel bandwidth (MHz) for OFDM modulation
        %   Specify the channel bandwidth for the packet as one of 'CBW20'
        %   | 'CBW10' | 'CBW5' to indicate 20MHz or 10MHz or 5MHz use
        %   respectively. The default is 'CBW20'.
        ChannelBandwidth = 'CBW20';     
                
        %MCS OFDM modulation and coding scheme
        %   Specify the modulation and coding scheme used for the
        %   transmission of a PPDU in the range of 0 to 7, inclusive, for
        %   OFDM modulation. The default is 0.
        MCS = 0;

        %DataRate Data rate in Mbps for DSSS modulation
        %   Specify the rate used to transmit the PSDU as one of the 
        %   following:
        %   '1Mbps'           - Differential Binary Phase Shift Keying
        %                       (DBPSK) modulation with 1Mbps data rate.
        %   '2Mbps'           - Differential Quadrature Phase Shift Keying
        %                       (DQPSK) modulation with 2Mbps data rate.
        %   '5.5Mbps'         - Complementary Code Keying (CCK) modulation
        %                       with 5.5Mbps data rate.
        %   '11Mbps'          - Complementary Code Keying (CCK) modulation
        %                       with 11Mbps data rate.
        % The default is '1Mbps'.
        DataRate = '1Mbps';
                
        %Preamble Preamble type for DSSS modulation
        %   Specify the PLCP preamble type as one of 'Long' | 'Short' to
        %   indicate use of the long PLCP preamble and short PLCP preamble
        %   respectively. For HR/DSSS modulation (Clause 17), the short
        %   PLCP preamble corresponds to the HR/DSSS/short mode. The
        %   default is 'Long'.
        Preamble = 'Long';
        
        %LockedClocks Clock locking indication for DSSS modulation
        %   Specify the "locked clocks bit" (b2) of the SERVICE field as
        %   one of true | false to indicate if the PHY implementation has
        %   its transmit frequency and symbol clocks derived from the same
        %   oscillator. For ERP-DSSS/CCK modulation, the PHY standard
        %   (Clause 19.1.3) requires that the implementation has locked
        %   clocks (b2 = 1) therefore LockedClocks should be set to true.
        %   The default is true.
        LockedClocks = true;
        
        %PSDULength PSDU length
        %   Specify the PSDU length, in bytes, for the data carried in a
        %   packet. The default is 1000 bytes.
        PSDULength = 1000;
        
        %NumTransmitAntennas Number of transmit antennas
        %   Specify the number of transmit antennas as an integer between 1
        %   and 8, inclusive.  The default is 1.
        NumTransmitAntennas = 1;
    end
    
    properties(Constant, Hidden)
        Modulation_Values = {'OFDM', 'DSSS'};
        ChannelBandwidth_Values  = {'CBW20', 'CBW10', 'CBW5'};
        DataRate_Values = {'1Mbps', '2Mbps', '5.5Mbps', '11Mbps'};
        Preamble_Values = {'Long', 'Short'};
    end
    
    methods
        function obj = wlanNonHTConfig(varargin)
            % Add sets for enum properties with different sized values
            obj = obj@wlanConfigBase('ChannelBandwidth', 'CBW20',  ...
                                     'DataRate', '1Mbps', ...
                                     'Preamble', 'Long', varargin{:});
        end
        
        % Self-validate and set properties
        function obj = set.Modulation(obj,val)
            prop = 'Modulation';
            validateEnumProperties(obj, prop, val);
            obj.(prop) = 'OFDM'; % default
            obj.(prop) = val;
        end
        
        function obj = set.ChannelBandwidth(obj,val)
            prop = 'ChannelBandwidth';
            validateEnumProperties(obj, prop, val);
            obj.(prop) = 'CBW20';    
            obj.(prop) = val;    
        end
        
        function obj = set.DataRate(obj,val)
            prop = 'DataRate';
            validateEnumProperties(obj, prop, val);
            obj.(prop) = '1Mbps'; % default
            obj.(prop) = val;
        end
        
        function obj = set.MCS(obj,val)
            prop = 'MCS';
            validateattributes(val, {'numeric'}, ...
                {'real','integer','scalar','>=',0,'<=',7}, ...
                [class(obj) '.' prop], prop);
            obj.(prop) = val;
        end
        
        function obj = set.Preamble(obj,val)
            prop = 'Preamble';
            validateEnumProperties(obj, prop, val);
            obj.(prop) = 'Long'; % default
            obj.(prop) = val;
        end
        
        function obj = set.LockedClocks(obj,val)
            prop = 'LockedClocks';
            validateattributes(val, {'logical'}, {'scalar'}, ...
                [class(obj) '.' prop], prop);
            obj.(prop)= val;
        end
        
        function obj = set.PSDULength(obj,val)
            prop = 'PSDULength';
            validateattributes(val, {'numeric'}, ...
                {'real','integer','scalar','>=',1,'<=',4095}, ...
                [class(obj) '.' prop], prop);
            obj.(prop) = val;
        end
        
        function obj = set.NumTransmitAntennas(obj, val)
            prop = 'NumTransmitAntennas';
            validateattributes(val, {'numeric'}, ...
                {'real','integer','scalar','>=',1,'<=',8}, ...
                [class(obj) '.' prop], prop);
            obj.(prop) = val;
        end
                                
        function varargout = validateConfig(obj, varargin)
            % validateConfig Validate the wlanNonHTConfig object.
            %
            %   validateConfig(CFGNONHT) validates the dependent properties
            %   for the specified wlanNonHTConfig parameter object.
            
            %   For INTERNAL use only, subject to future changes:
            %
            %   validateConfig(CFGNONHT, MODE) validates only the subset of 
            %   dependent properties as specified by the MODE input.
            %   MODE must be one of: TBD
            
            narginchk(1,2);nargoutchk(0,1);
            if (nargin==2)
                mode = varargin{1};
            else
                mode = 'Full';
            end
            
            s = struct('NumDataSymbols', -1, ...
                'NumPadBits',     -1, ...
                'NumPPDUSamples', -1, ...
                'TxTime',         -1);

            if strcmp(obj.Modulation, 'OFDM')
                if strcmp(mode, 'Full')
                    s = validate(obj);
                end
            % else
                % None needed for DSSS
            end
            
            if nargout==1
                varargout{1} = s;
            end
            
        end % end of validation function
                
    end
    
    methods (Access=protected)
        function flag = isInactiveProperty(obj, prop)
            % Controls the conditional display of properties

            flag = false;

            % ChannelBandwidth only for OFDM
            if strcmp(prop, 'ChannelBandwidth')
                flag = ~strcmp(obj.Modulation, 'OFDM');
            end

            % NumTransmitAntennas only for OFDM and for CBW20
            if strcmp(prop, 'NumTransmitAntennas')
                flag = ~strcmp(obj.Modulation, 'OFDM') || ...
                       ~strcmp(obj.ChannelBandwidth, 'CBW20');
            end
            
            % MCS only for OFDM
            if strcmp(prop, 'MCS')
                flag = ~strcmp(obj.Modulation, 'OFDM');
            end
            
            % DataRate only for DSSS            
            if strcmp(prop, 'DataRate')
                flag = ~strcmp(obj.Modulation, 'DSSS');
            end
            
            % Preamble only for DSSS and for DataRate > 1 Mbps
            if strcmp(prop, 'Preamble')
                flag = ( ~strcmp(obj.Modulation, 'DSSS') || ...
                         ( strcmp(obj.Modulation, 'DSSS') && ...
                           strcmp(obj.DataRate, '1Mbps') ) );
            end
            
            % LockedClocks only for DSSS
            if strcmp(prop, 'LockedClocks')
                flag = ~strcmp(obj.Modulation, 'DSSS');
            end
            
        end
    end
    
    methods (Access=private)
        function s = privInfo(obj)
            %privInfo Returns information relevant to the object
            %   S = privInfo(CFGNONHT) returns a structure, S, containing
            %   the relevant information for the wlanNonHTConfig object,
            %   CFGNONHT. Only OFDM modulation type is supported.
            %
            %   The output structure S has the following fields:
            %
            %   NumDataSymbols - Number of OFDM symbols for the Data field
            %   NumPadBits     - Number of pad bits in the Data field
            %   NumPPDUSamples - Number of PPDU samples per transmit antenna
            %   TxTime         - The time in microseconds, required to
            %                    transmit the PPDU.
            
            s = struct('NumDataSymbols', -1, ...
                'NumPadBits',     -1, ...
                'NumPPDUSamples', -1, ...
                'TxTime',         -1);
            
            % Compute the number of OFDM symbols in Data field
            mcsTable  = getRateTable(obj);
            numDBPS = mcsTable.NDBPS;
            
            Ntail = 6; Nservice = 16;
            numDataSym = ceil((8*obj.PSDULength + Nservice + Ntail)/numDBPS);
            numPadBits = numDataSym * numDBPS - (8*obj.PSDULength + Nservice + Ntail);
            
            % Compute the number of PPDU samples at CBW
            Nfft = 64;      % Fixed
            cpLen = Nfft/4; % Always long
            
            numSymbols = 2 + 2 + 1 + numDataSym;
            numSamples = numSymbols*(Nfft + cpLen);
            switch obj.ChannelBandwidth
                case 'CBW10'
                    txTime = numSymbols*8;
                case 'CBW5'
                    txTime = numSymbols*16;
                otherwise % CBW20
                    txTime = numSymbols*4;
            end
            
            % Set output structure
            s.NumDataSymbols = numDataSym;
            s.NumPadBits     = numPadBits;
            s.NumPPDUSamples = numSamples;
            s.TxTime         = txTime;            
        end
        
        function s = validate(obj)
            s = privInfo(obj);
        end
            
    end
end

% [EOF]