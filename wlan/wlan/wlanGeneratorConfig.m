classdef wlanGeneratorConfig < wlanConfigBase 
%   Warning: The use of wlanGeneratorConfig object is discouraged for
%   parameterizing the <a
%   href="matlab:help('wlanWaveformGenerator')">wlanWaveformGenerator</a> function. See the 
%   documentation of wlanWaveformGenerator for the recommended parameter
%   name-value pair syntax.
%   
%wlanGeneratorConfig Waveform generator configuration
%   CFGWAVEGEN = wlanGeneratorConfig creates a waveform generator
%   configuration object, CFGWAVEGEN. This object is used to configure 
%   <a href="matlab:help('wlanWaveformGenerator')">wlanWaveformGenerator</a>.
%
%   CFGWAVEGEN = wlanGeneratorConfig(Name,Value) creates a waveform
%   generator configuration object, CFGWAVEGEN, with the specified property
%   NAME set to the specified VALUE. Additional name-value pair arguments
%   can be specified in any order as (Name1,Value1,...,NameN,ValueN).
%
%   wlanGeneratorConfig properties:
%   NumPackets              - Number of packets to generate.
%   IdleTime                - Time in seconds after each packet.
%   ScramblerInitialization - The initial state of the data scrambler for
%                             each packet to generate. This is only
%                             applicable for OFDM based formats.
%   Windowing               - Enable or disable windowing. This is only 
%                             applicable for OFDM based formats.
%   WindowTransitionTime    - Window length in seconds, applicable when
%                             windowing is enabled.
%
%   Example:
%   %  Configure wlanWaveformGenerator to produce a time domain signal
%   %  txWaveform for an 802.11ac VHT transmission with 10 packets and 20
%   %  microsecond idle period between packets. A random scrambler initial
%   %  value is used for each packet.
%   
%   cfgVHT = wlanVHTConfig;             % Create format configuration   
%   
%   numPackets    = 10;                 % Generate 10 packets
%   idleTime      = 20e-6;              % 20 microsecond idle period
%   scramblerInit = randi([1 127], ...  % Random scrambler initial state
%                   numPackets, 1); 
%         
%   % Method 1: Use wlanGeneratorConfig object in wlanWaveformGenerator
%   cfgWaveGen = wlanGeneratorConfig(....
%       'NumPackets',              numPackets, ...
%       'IdleTime',                idleTime, ...
%       'ScramblerInitialization', scramblerInit);
%
%   txWaveform1 = wlanWaveformGenerator([1;0;0;1], cfgVHT, cfgWaveGen);
%           
%   % Mehothod 2: Use name-value pair syntax in wlanWaveformGenerator
%   txWaveform2 = wlanWaveformGenerator([1;0;0;1], cfgVHT, ...
%       'NumPackets',              numPackets, ...
%       'IdleTime',                idleTime, ...
%       'ScramblerInitialization', scramblerInit);
%
%   See also wlanWaveformGenerator.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

% Public properties
properties (Access = public)
    %NumPackets Number of packets
    %   Specify the number of packets to generate. The default is 1.
    NumPackets = 1;

    %IdleTime Idle time after each packet
    %   Specify the length in seconds of an idle period included between
    %   multiple generated packets and after the last packet generated.
    %   IdleTime must be 0 or greater than 2e-6 seconds. The default is 0
    %   seconds.
    IdleTime = 0; % 0 seconds

    %ScramblerInitialization Scrambler initial states
    %   Specify the scrambler initial states for OFDM based formats as a
    %   double or int8-typed scalar, or an Np-by-Nu matrix, containing
    %   integer values between 1 and 127 inclusive. Np is the number of
    %   packets to generate and Nu is the number of users. If a scalar is
    %   provided all packets are initialized with the same state for all
    %   users. Specifying a matrix allows a different initial state to be
    %   used per user and per packet. Each column specifies the initial
    %   states for a single user, therefore up to 4 columns are supported.
    %   If a single column is provided, the same initial states will be
    %   used for all users. Each row represents the initial state of each
    %   packet to generate. Therefore a matrix with multiple rows allows a
    %   different initial state to be used per packet, where the first row
    %   contains the initial state of the first packet. Internally the rows
    %   are looped if the number of packets to generate exceeds the number
    %   of rows of the matrix provided. The default value of this property
    %   is 93 which is the example state given in IEEE Std 802.11-2012
    %   Section L.1.5.2.
    ScramblerInitialization = 93;

    %Windowing Enable windowing
    %   Set this property to true to enable windowing between consecutive
    %   OFDM symbols. The default is true.
    Windowing = true;

    %WindowTransitionTime Windowing transition time
    %   Specify the length of the windowing transition in seconds applied
    %   to OFDM symbols when windowing is enabled. WindowTransitionTime is
    %   a positive scalar with a maximum value of 6.4e-06 seconds. No
    %   windowing is applied for a transition time of 0. The default is
    %   1.0e-07 seconds.
    WindowTransitionTime = 1.0e-07;
end

methods 
    % Constructor allows setting properties with name-value pairs        
    function obj = wlanGeneratorConfig(varargin)
      coder.internal.warning('wlan:wlanGeneratorConfig:Deprecation');
      obj = obj@wlanConfigBase(varargin{:});
    end

    % Validate NumPackets is a scalar numeric integer >= 0
    function obj = set.NumPackets(obj,val)
        validateattributes(val,{'numeric'},{'scalar','integer','>=',0}, ...
            'wlanGeneratorConfig','NumPackets');
        obj.NumPackets = val;
    end  

    % Validate IdleTime is a scalar numeric >= 0
    function obj = set.IdleTime(obj,val)
        validateattributes(val,{'numeric'}, ...
            {'scalar','real','>=',0},'wlanGeneratorConfig','IdleTime');

        coder.internal.errorIf(val > 0 && val < 2e-6 , ...
        'wlan:wlanConfig:InvalidIdleTimeValue' );
        obj.IdleTime = val;

    end 

    % Validate ScramblerInitialization is a scalar numeric integer >= 1 <
    % 2^7 (7 bits)
    function obj = set.ScramblerInitialization(obj, val)
        validateattributes(val,{'double','int8'}, ...
            {'real','integer','2d','nonempty','>=',1,'<=',127}, ...
            'wlanGeneratorConfig','ScramblerInitialization');

        coder.internal.errorIf(size(val, 2) > 4, ...
            'wlan:wlanGeneratorConfig:ScramInitMoreThan4Col');
        obj.ScramblerInitialization = val;
    end 

    % Validate Windowing 
    function obj = set.Windowing(obj,val)
        validateattributes(val, {'logical'}, {'scalar'}, ...
            'wlanGeneratorConfig','Windowing');
        obj.Windowing= val;
    end

    % Validate WindowTransitionTime 
    function obj = set.WindowTransitionTime(obj,val)
        validateattributes(val, {'numeric'}, ...
            {'real','scalar','>=',0,'<=',(6.4e-6)}, ...
             'wlanGeneratorConfig','WindowTransitionTime');
        obj.WindowTransitionTime = val;
    end
end               

methods (Access = protected)
    function flag = isInactiveProperty(obj, prop)
        flag = false;
        if strcmp(prop, 'WindowTransitionTime')
            flag = ~(obj.Windowing);
        end
    end
end
end