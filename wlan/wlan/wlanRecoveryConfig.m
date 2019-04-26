classdef wlanRecoveryConfig < wlanConfigBase    
%wlanRecoveryConfig Construct a configuration object for data recovery
%   CFGREC = wlanRecoveryConfig constructs a configuration object for
%   recovering the data in signaling or data fields. Adjust the property
%   values of the object, which indicate different algorithm parameters or
%   operations at the receiver, to achieve optimal recovery performance.
%
%   CFGREC = wlanRecoveryConfig(Name,Value) constructs a recovery
%   configuration object, CFGREC, with the specified property Name set to
%   the specified Value. You can specify additional name-value pair
%   arguments in any order as (Name1,Value1,...,NameN,ValueN).
% 
%   wlanRecoveryConfig properties:
%
%   OFDMSymbolOffset   - OFDM symbol sampling offset
%   EqualizationMethod - Equalization method
%   PilotPhaseTracking - Pilot phase tracking
% 
%   % Example: 
%   %    Create a wlanRecoveryConfig object for performing ZF equalization 
%   %    and OFDM symbol sampling offset of 0.5 in a recovery process
% 
%   cfgRec = wlanRecoveryConfig( ...
%       'OFDMSymbolOffset',   0.5, ...
%       'EqualizationMethod', 'ZF')
%  
%   See also wlanLSIGRecover, wlanVHTSIGARecover, wlanVHTSIGBRecover,
%   wlanVHTDataRecover, wlanHTSIGRecover, wlanHTDataRecover,
%   wlanNonHTDataRecover.

% Copyright 2015 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

properties (Access = 'public')
    %OFDMSymbolOffset OFDM symbol offset
    %   Specify the sampling offset as a fraction of the cyclic prefix (CP)
    %   length for every OFDM symbol, as a double precision, real scalar
    %   between 0 and 1, inclusive. The OFDM demodulation is performed
    %   based on Nfft samples following the offset position, where Nfft
    %   denotes the FFT length. The default value of this property is 0.75,
    %   which means the offset is three quarters of the CP length.
    OFDMSymbolOffset = 0.75;
    %EqualizationMethod Equalization method
    %   Specify the equalization method as one of 'MMSE' | 'ZF'. The
    %   default value of this property is 'MMSE'.
    EqualizationMethod = 'MMSE';
    %PilotPhaseTracking Pilot phase tracking
    %   Specify the pilot phase tracking performed as one of 'PreEQ' |
    %   'None'. 'PreEQ' pilot phase tracking estimates and corrects a
    %   common phase offset across all subcarriers and receive antennas for
    %   each received OFDM symbol before equalization. The default is
    %   'PreEQ'.
    PilotPhaseTracking = 'PreEQ';
end

properties(Constant, Hidden)
    EqualizationMethod_Values = {'MMSE', 'ZF'};
    PilotPhaseTracking_Values = {'PreEQ', 'None'};
end

methods
  function obj = wlanRecoveryConfig(varargin)
    obj = obj@wlanConfigBase('EqualizationMethod', 'MMSE', ...
        'PilotPhaseTracking','PreEQ', varargin{:});
  end
  
  function obj = set.OFDMSymbolOffset(obj, val)
    prop = 'OFDMSymbolOffset';
    validateattributes(val, {'double'}, ...
        {'real','scalar','>=',0,'<=',1}, ...
        [class(obj) '.' prop], prop);
    obj.(prop) = val;
  end
  
  function obj = set.EqualizationMethod(obj, val)
    prop = 'EqualizationMethod';
    validateEnumProperties(obj, prop, val);
    obj.(prop) = ''; 
    obj.(prop) = val;
  end
  
  function obj = set.PilotPhaseTracking(obj,val)
    prop = 'PilotPhaseTracking';
    validateEnumProperties(obj, prop, val);
    obj.(prop) = ''; 
    obj.(prop) = val;
  end
end

end

% [EOF]