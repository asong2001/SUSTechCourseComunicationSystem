classdef (Hidden) ChannelBase < matlab.System & matlab.system.mixin.Propagates
%ChannelBase Filter input signal through a TGn and TGac fading channel
%   The implementation is based on TGn and TGac channel model, as specified
%   by the IEEE 802.11 Wireless LAN Working group [1,2]. The power delay
%   profile, spatial correlation and Doppler filter coefficients values are
%   based on TGn Implementation note version 3.2 - May 2004.
%
%   %   References:
%   [1] Erceg, V. et al. "TGn Channel Models."  Doc. IEEE802.11-03/940r4.
%   [2] Briet, G. et al. "TGac Channel Model Addendum." Doc. IEEE 802.11-09/0308r12

% Copyright 2015 The MathWorks, Inc.

%#codegen
%#ok<*EMCLS>

properties (Nontunable)
    %DelayProfile Channel propagation condition
    %   Specify the Delay profile model of WLAN multipath fading channel as
    %   one of 'Model-A' | 'Model-B' | 'Model-C' | 'Model-D' | 'Model-E' |
    %   'Model-F'. The default is 'Model-B'.
    DelayProfile = 'Model-B';
    % CarrierFrequency Carrier frequency (Hz) 
    %   Specify the carrier frequency of the input signal in Hz. The
    %   default is 5.25e9 Hz.
    CarrierFrequency = 5.25e9;
    %LargeScaleFadingEffect Large scale fading effect
    %   Specify the large scale fading effect in WLAN multipath fading
    %   channel as one of 'None' | 'Pathloss' | 'Shadowing' | 'Pathloss and
    %   shadowing'. The default is 'None'.
    LargeScaleFadingEffect = 'None';
    %TransmitReceiveDistance Distance between transmitter and receiver (m)
    %   Specify the separation between transmitter and receiver in meters.
    %   Used to compute the path loss, and to determine whether the channel
    %   has LOS or NLOS condition. The path loss and standard deviation of
    %   shadow fading loss is dependent on the separation between the
    %   transmitter and the receiver. The default is 3.
    TransmitReceiveDistance = 3;
    %PowerLineFrequency Power line frequency (Hz)
    %   Specify the power line frequency as one of '50Hz' | '60Hz'. The
    %   power line frequency is 60Hz in US and 50Hz in Europe. This
    %   property only applies when DelayProfile is set to Model-D or
    %   Model-E. The default is 60Hz.
    PowerLineFrequency = '60Hz';
    %RandomStream Source of random number stream
    %   Specify the source of random number stream as one of 'Global
    %   stream' | 'mt19937ar with seed'. If RandomStream is set to 'Global
    %   stream', the current global random number stream is used for
    %   normally distributed random number generation, in which case the
    %   reset method only resets the filters. If RandomStream is set to
    %   'mt19937ar with seed', the mt19937ar algorithm is used for normally
    %   distributed random number generation, in which case the reset
    %   method not only resets the filters but also re-initializes the
    %   random number stream to the value of the Seed property. The default
    %   value of this property is 'Global stream'.
    RandomStream = 'Global stream';  
    %Seed Initial seed
    %   Specify the initial seed of a mt19937ar random number generator
    %   algorithm as a double precision, real, nonnegative integer scalar.
    %   This property applies when you set the RandomStream property to
    %   'mt19937ar with seed'. The Seed is to re-initialize the mt19937ar
    %   random number stream in the reset method. The default value of this
    %   property is 73.
    Seed = 73;
end

properties (Nontunable, Logical)
    %NormalizePathGains Normalize average path gains to 0 dB
    %   Set this property to true to normalize the fading processes such
    %   that the total power of the path gains, averaged over time, is 0dB.
    %   The default is true.
    NormalizePathGains = true;
    %FluorescentEffect Enable fluorescent effect
    %   Set this property to add fluorescent effect. This property is only
    %   applicable for Delayprofile, Model-D and Model-E. The default is
    %   true.
    FluorescentEffect = true;
    %NormalizeChannelOutputs Normalize output by number of receive antennas
    %   Set this property to true to normalize the channel outputs by the
    %   number of receive antennas. The default is true.
    NormalizeChannelOutputs = true;
end

properties (Access = private, Nontunable)
    % Number of transmit antennas
    pNt;
    % Number of receive antennas
    pNr;
    % Number of links, NL = pNt * pNr
    pNumLinks;
    % Channel sampling frequency
    pFc;
    % Normalized channel sampling rate
    pDoppler; 
    % Doppler filter weights
    pNumeratorB;
    pDenominatorB;
    pNumeratorS;
    pDenominatorS;
end

properties (Access = private)
    % Number of multipaths
    pNumPaths;
    % Rician K factor in decibels
    pKfactor;
    % Path power
    pPathPower;
    % Antenna correlation
    pAntennaCorrelation;
    % Rician LOS component
    pRicianLOS;
    % Rician NLOS component
    pRicianNLOS;
    % Channel filter tap gain sinc interpolation matrix
    pSincInterpMatrix;
    % Number of channel taps
    pNumChannelTaps;
    % Running Time Instance
    pRunningTimeInstant;
    % Index for oversampled channel samples
    pOversampledFadingTime;
    % Matrix to store the fadings samples before interpolation.
    pOversampledFadingMatrix;
    % Vector of channel samples at channel sampling rate
    pFadingTime;
end

properties (Access = private)
    % Seeding state
    pWGNState;
    % Filter state
    pFilterState;
    % Random phase
    pRandomInterfererToCarrierRatio;
    % Filter state
    pChannelFilterState;
    % State retention
    %pLastTwoColumn
    pLastTwoColumn;
    % Fading offset
    pFadingOffset;
    % Looping variables
    pRunningTime;
    % Fading sample reserve index
    pFadingBufferIndex;
    % Fading Buffer
    pFadingBuffer;
    % Dummy offset
    pOffset;
    % Fluorescent effect random phase
    pRandomPhase;
    % Large scale losses
    pLargeScaleLosses;
    % Stream
    pStream;
    % Power line frequency
    pPowerLineFrequency;
    % Default shadowing value
    pDefaultShadowFading;
end

% Child classes TGnChannel and TGacChannel can access these properties
properties (Access = protected)
    % Channel filter delay, measured in samples
    pChannelFilterDelay = 0;
    % Path power
    pPathPowerdBs;
    % Path delays
    pPathDelays;
    % Legacy generator
    pLegacyGenerator;
    % ShadowFading in decibels
    pShadowFading = 0;
    % Path loss in decibels
    pPathloss = 0;
end

properties(Constant, Hidden)
    RandomStreamSet = matlab.system.StringSet({'Global stream', 'mt19937ar with seed'});
    DelayProfileSet = matlab.system.StringSet({'Model-A','Model-B','Model-C','Model-D','Model-E','Model-F'});
    LargeScaleFadingEffectSet= matlab.system.StringSet({'None','Pathloss','Shadowing','Pathloss and shadowing'});
    PowerLineFrequencySet = matlab.system.StringSet({'50Hz','60Hz'});
end

methods
  function obj = ChannelBase(varargin) % Constructor
    setProperties(obj, nargin, varargin{:});
    obj.pPathDelays = 0;
    obj.pPathPowerdBs = 0;
  end 
    
   function set.CarrierFrequency(obj,val)
     propName = 'CarrierFrequency';
        validateattributes(val, {'numeric'}, ...
            {'real','scalar','>',0,'finite'}, ...
            [class(obj) '.' propName], propName); 
      obj.CarrierFrequency= val;
   end
   
  function set.TransmitReceiveDistance(obj, val)
    propName = 'TransmitReceiveDistance';
    validateattributes(val, {'numeric'}, ...
        {'real','scalar','>',0,'finite'}, ...
        [class(obj) '.' propName], propName);    
    obj.TransmitReceiveDistance = val;
  end
    
  function set.Seed(obj, seed)
    propName = 'Seed';
    validateattributes(seed, {'double'}, ...
        {'real','scalar','integer','nonnegative','finite'}, ...
        [class(obj) '.' propName], propName);  %#ok<*EMCA
    obj.Seed = seed;
  end 
end

methods(Access = protected)
  function num = getNumInputsImpl(~)
    num = 1;
  end
  
  function num = getNumOutputsImpl(~)
    num = 2;
  end  
  
  function validateInputsImpl(obj, x)    
      validateattributes(x, {'double'}, ... 
          {'2d', 'finite', 'ncols', obj.NumTransmitAntennas}, ...
          [class(obj) '.' 'Signal input'], 'Signal input');     
   end
     
  function setupImpl(obj, varargin)
     
    % Initialization of TGn and TGac channel model parameters
    getInfoParameters(obj);
    
    % Set power line frequency
    if strcmp(obj.PowerLineFrequency,'50Hz')
        obj.pPowerLineFrequency = 50;
    else
        obj.pPowerLineFrequency = 60;
    end
    
    coder.internal.errorIf((10*obj.pPowerLineFrequency)>(5*obj.pFc), ...
        'wlan:wlanChannel:AliasingfluorescentFrequency');
    
    
    % Computation of PDP for LOS and NLOS scenario
    [obj.pRicianLOS, obj.pRicianNLOS, obj.pPathPower] ...
                                    = ricianComponents(obj);
    
    % Setup filter weights (Doppler filter coefficients)
    obj.pNumeratorB = [ 2.785150513156437e-04
                  -1.289546865642764e-03
                   2.616769929393532e-03
                  -3.041340177530218e-03
                   2.204942394725852e-03
                  -9.996063557790929e-04
                   2.558709319878001e-04
                  -2.518824257145505e-05];
    obj.pDenominatorB =[1.000000000000000
                  -5.945307133332568
                   1.481117656568614e+01
                  -1.985278212976179e+01
                   1.520727030904915e+01
                  -6.437156952794267
                   1.279595585941577
                  -6.279622049460144e-02];

    obj.pNumeratorS= [7.346653645411014e-02
    -2.692961866120396e-01-3.573876752780644e-02i
     5.135568284321728e-01+1.080031147197969e-01i
    -5.941696436704196e-01-1.808579163219880e-01i
     4.551285303843579e-01+1.740226068499579e-01i
    -2.280789945309740e-01-1.111510216362428e-01i
     7.305274345216785e-02+4.063385976397657e-02i
    -1.281829045132903e-02-9.874529878635624e-03i];

    obj.pDenominatorS=[1.0
    -4.803186161814090e+00-6.564741337993197e-01i
     1.083699602445439e+01+2.655291948209460e+00i
    -1.461669273797422e+01-5.099798708358494e+00i
     1.254720870662511e+01+5.723905805435745e+00i
    -6.793962089799931e+00-3.887940726966912e+00i
     2.107762007599925e+00+1.503178356895452e+00i
    -2.777250532557400e-01-2.392692168502482e-01i];
   
  end
  
  function resetImpl(obj)
   
    if coder.target('MATLAB') 
       if strcmp(obj.RandomStream, 'Global stream')
           obj.pStream = RandStream.getGlobalStream;
       else
           obj.pStream = RandStream('mt19937ar', 'seed', obj.Seed);
           obj.pWGNState = obj.pStream.State;
       end
    elseif strcmp(obj.RandomStream, 'mt19937ar with seed')
          obj.pStream = eml_rand_mt19937ar('preallocate_state');
          obj.pWGNState= eml_rand_mt19937ar('seed_to_state', ...
                           obj.pStream, obj.Seed);
    end
       
  
   wgnoise = generateWGNoise(obj);
   obj.pShadowFading = double(wgnoise(1,1).*obj.pDefaultShadowFading);
      
   % Initialize Doppler filter
   filterOrder = 7;
   obj.pFilterState = complex(zeros(filterOrder, ...
                                    obj.pNt*obj.pNr*obj.pNumPaths));
                                
   % Initialize the filter in order to avoid transient state
   temp = generateNoiseSamples(obj, 1000); %#ok<NASGU>
     
   % Initialization of fluorescent effects random variables
   GaussianStandardDeviation = 0.0107;
   GaussianMean              = 0.0203;

   % Gaussian noise generation
   wgnoise = generateWGNoise(obj);

   obj.pRandomInterfererToCarrierRatio = ...,
        (double(wgnoise(1,1)*GaussianStandardDeviation + GaussianMean)).^2;

   % Set random phase
   if coder.target('MATLAB')   
     if strcmp(obj.RandomStream, 'Global stream') 
        stream = RandStream.getGlobalStream;
        wgnoise = rand(stream,1,1);
     else 
        stream = RandStream('mt19937ar', 'seed', obj.Seed);
        wgnoise = rand(stream,1,1);
     end
   else
     if strcmp(obj.RandomStream, 'Global stream') 
         wgnoise = rand(1,1);
     else
        stream = eml_rand_mt19937ar('preallocate_state');
        state = eml_rand_mt19937ar('seed_to_state', stream, obj.Seed);
        for idx = 1:1
           [state, wgnoise] = eml_rand_mt19937ar('generate_uniform',state);
        end
     end 
   end
   
   obj.pRandomPhase = double(wgnoise(1,1))*2*pi;
   
   if obj.pNumChannelTaps > 1
      obj.pChannelFilterState = complex(zeros(obj.pNumChannelTaps-1, ...
                                obj.pNt*obj.pNumPaths));
   else
      obj.pChannelFilterState = complex(0);
   end
   
   % Setup large scale fading effect for the simulation
    if ~strcmp(obj.LargeScaleFadingEffect,'None')
       switch obj.LargeScaleFadingEffect
         case 'Shadowing'
            obj.pLargeScaleLosses = 1/sqrt(10.^(0.1.*(obj.pShadowFading)));
         case 'Pathloss'
            obj.pLargeScaleLosses = 1/sqrt(10.^(0.1.*(obj.pPathloss)));
         otherwise
            obj.pLargeScaleLosses = 1/sqrt(10.^(0.1.*(obj.pPathloss ...
                                                + obj.pShadowFading)));
       end
    end
    
   % Initialization 
   obj.pOffset            = 0;
   obj.pLastTwoColumn     = complex(zeros( ...
                                   [obj.pNumPaths*obj.pNumLinks, 2]));
   obj.pRunningTimeInstant = 0;
   obj.pFadingOffset = 0;
   obj.pOversampledFadingTime = zeros(1,1);
   obj.pOversampledFadingMatrix = complex(zeros(1,1));
   obj.pFadingTime = zeros(1,1);

  end
  
  function varargout = stepImpl(obj, input)
   
   Ns = size(input,1);
   nTx = obj.pNt;
   nRx = obj.pNr;
   numPaths = obj.pNumPaths;

   % Output empty signals and path gains for an empty input
   if Ns == 0
        varargout{1} = zeros(Ns, nRx);
        varargout{2} = NaN(Ns, numPaths, nTx, nRx);
        return;
   end
   
   % Number of samples in time
   chSampleTime = 1/obj.pFc;
   inpSampleTime = 1/obj.SampleRate;
   Nt = Ns* inpSampleTime;
              
   runningTime = obj.pRunningTimeInstant:inpSampleTime: ...
                    (Ns-1)*inpSampleTime + obj.pRunningTimeInstant;
   
   if runningTime(1)==0
       numChSamples = ceil(Nt/chSampleTime) + 1;                 
       pathGains = generateFadingSamples(obj,runningTime,numChSamples);
       obj.pOffset = 1;
   else
       if runningTime(end) < obj.pOversampledFadingTime(end)
           % No additional channel samples are required to process the input samples.
           pathGains = (interp1(obj.pOversampledFadingTime.', ...
             obj.pOversampledFadingMatrix.',runningTime.','linear'));
       
       else
           numInterpolatedSamples = size(runningTime,2);
           % Use the remaining samples before generating more.
           previousIdx = find(runningTime<obj.pFadingTime(end));

           % Generate output samples before generating new channel samples
           H1 = (interp1(obj.pOversampledFadingTime.', ...
                obj.pOversampledFadingMatrix.', ...
                runningTime(previousIdx).','linear'));
           
           % Additional samples needed to process the input.
           newInterpolatedSamples = numInterpolatedSamples - size(previousIdx,2);
        
           % Time required to process the remaining interpolated samples.
           timeInterpolatedSamples = newInterpolatedSamples * inpSampleTime;
           
           % Number of channel samples required to process this input
           numChannelSamples = ceil(timeInterpolatedSamples/chSampleTime) + 1;

           % The isempty check is to avoid the special case when number of
           % interpolated channel samples are exact match to the length of
           % samples between two channel samples
           runningTimeInstant = runningTime(1,size(runningTime,2) ...
                                + isempty(previousIdx) ...
                                -newInterpolatedSamples) + ...
                                inpSampleTime*(1-isempty(previousIdx));
                   
           runningTime = runningTimeInstant:inpSampleTime: ...
               (newInterpolatedSamples-1)*inpSampleTime + runningTimeInstant;
          
           H2 = generateFadingSamples(obj,runningTime,numChannelSamples);
           pathGains = [H1;H2];
       end
   end
       
   % Update states
   obj.pFadingOffset = obj.pFadingTime(size(obj.pFadingTime,2))+ ...
                           chSampleTime;
   obj.pRunningTimeInstant  = runningTime(1,size(runningTime,2)) + inpSampleTime;
       
   varargout{1} = channelFilter(obj, input, pathGains);
 
   if obj.NormalizeChannelOutputs
         %Normalize by the number of selected receive antennas so that the
         %total output power is equal to the total input power
         varargout{1} = varargout{1}/sqrt(nRx);
   end

   varargout{2} = permute(reshape(pathGains, Ns, nRx, nTx, ...
                                              numPaths), [1, 4, 3, 2]);
    
  end
  
  function flag = isInactivePropertyImpl(obj, prop)
    if strcmp(prop, 'FluorescentEffect') ||  strcmp(prop, 'PowerLineFrequency') 
        flag = ~(strcmpi(obj.DelayProfile, 'Model-D') || strcmpi(obj.DelayProfile, 'Model-E'));
    elseif strcmp(prop, 'Seed')
        flag = strcmp(obj.RandomStream, 'Global stream');
    else 
        flag = false;
    end
  end
  
  function releaseImpl(~)

  end
  
  function s = saveObjectImpl(obj)
    s = saveObjectImpl@matlab.System(obj);
    if isLocked(obj)
        s.pNt = obj.pNt;
        s.pNr = obj.pNr;
        s.pNumPaths = obj.pNumPaths; 
        s.pPathPower = obj.pPathPower;
        s.pPathPowerdBs = obj.pPathPowerdBs;
        s.pNumLinks = obj.pNumLinks;
        s.pNumChannelTaps = obj.pNumChannelTaps; 
        s.pSincInterpMatrix = obj.pSincInterpMatrix;
        s.pWGNState = obj.pWGNState;
        s.pChannelFilterState = obj.pChannelFilterState; 
        s.pRunningTime = obj.pRunningTime; 
        s.pOffset = obj.pOffset;
        s.pFc = obj.pFc;
        s.pRandomPhase = obj.pRandomPhase;
        s.pRicianLOS = obj.pRicianLOS;
        s.pRicianNLOS = obj.pRicianNLOS;
        s.pPathPower = obj.pPathPower;
        s.pAntennaCorrelation = obj.pAntennaCorrelation;
        s.pPathloss = obj.pPathloss;
        s.pDoppler = obj.pDoppler;
        s.pShadowFading = obj.pShadowFading;
        s.pFilterState = obj.pFilterState;
        s.pLastTwoColumn = obj.pLastTwoColumn;
        s.pRunningTime = obj.pRunningTime;
        s.pRandomPhase = obj.pRandomPhase;
        s.pRandomInterfererToCarrierRatio = obj.pRandomInterfererToCarrierRatio;
        s.pLargeScaleLosses = obj.pLargeScaleLosses;
        s.pPowerLineFrequency = obj.pPowerLineFrequency;
        s.pNumeratorB = obj.pNumeratorB;
        s.pDenominatorB = obj.pDenominatorB;
        s.pNumeratorS = obj.pNumeratorS;
        s.pDenominatorS = obj.pDenominatorS;
        s.pLegacyGenerator = obj.pLegacyGenerator;
        s.pRunningTimeInstant = obj.pRunningTimeInstant;
        s.pOversampledFadingTime = obj.pOversampledFadingTime;
        s.pOversampledFadingMatrix = obj.pOversampledFadingMatrix;
        s.pFadingTime = obj.pFadingTime;
        s.pFadingOffset = obj.pFadingOffset;
        s.pChannelFilterDelay = obj.pChannelFilterDelay;
        s.pPathDelays = obj.pPathDelays;
        s.pDefaultShadowFading = obj.pDefaultShadowFading;
    end
  end
  
  function loadObjectImpl(obj, s, wasLocked)
    if wasLocked
        obj.pNt = s.pNt;
        obj.pNr = s.pNr;
        obj.pNumPaths = s.pNumPaths;
        obj.pPathPower = s.pPathPower;
        obj.pPathPowerdBs = s.pPathPowerdBs;
        obj.pNumLinks = s.pNumLinks;
        obj.pNumChannelTaps = s.pNumChannelTaps; 
        obj.pSincInterpMatrix = s.pSincInterpMatrix;
        obj.pWGNState = s.pWGNState;
        obj.pChannelFilterState = s.pChannelFilterState;
        obj.pRunningTime = s.pRunningTime;
        obj.pOffset = s.pOffset;
        obj.pFc = s.pFc;
        obj.pRandomPhase = s.pRandomPhase;
        obj.pRicianLOS = s.pRicianLOS;
        obj.pRicianNLOS = s.pRicianNLOS;
        obj.pPathPower = s.pPathPower;
        obj.pDoppler = s.pDoppler;
        obj.pAntennaCorrelation = s.pAntennaCorrelation;
        obj.pPathloss = s.pPathloss;
        obj.pShadowFading = s.pShadowFading;
        obj.pFilterState = s.pFilterState;
        obj.pLastTwoColumn = s.pLastTwoColumn;
        obj.pRunningTime = s.pRunningTime;
        obj.pRandomPhase = obj.pRandomPhase;
        obj.pRandomInterfererToCarrierRatio = s.pRandomInterfererToCarrierRatio;
        obj.pLargeScaleLosses = s.pLargeScaleLosses;
        obj.pPowerLineFrequency = s.pPowerLineFrequency;
        obj.pNumeratorB = s.pNumeratorB;
        obj.pDenominatorB = s.pDenominatorB;
        obj.pNumeratorS = s.pNumeratorS;
        obj.pDenominatorS = s.pDenominatorS;
        obj.pLegacyGenerator = s.pLegacyGenerator;
        obj.pRunningTimeInstant = s.pRunningTimeInstant;
        obj.pOversampledFadingTime = s.pOversampledFadingTime;
        obj.pOversampledFadingMatrix = s.pOversampledFadingMatrix;
        obj.pFadingTime = s.pFadingTime;
        obj.pFadingOffset = s.pFadingOffset;
        obj.pChannelFilterDelay = s.pChannelFilterDelay;
        obj.pPathDelays = s.pPathDelays;
        obj.pDefaultShadowFading = s.pDefaultShadowFading;
    end
    % Do not call the loadObjectImpl method of the matlab.System object
    % because it does not support private properties of the base class.
    % For example, pNumPaths is not a property of wlanTGnChannel and will
    % not be excluded from s before calling set(obj, s).
    field = fieldnames(s);
    for i = 1 : length(field)
        curField = field{i};
        if strcmp(curField(1), 'p')
            s = rmfield(s, curField);
        end
    end
    set(obj, s);
  end
  
  function flag = isInputSizeLockedImpl(~,~)
     flag = false;
  end
  
  function varargout = isOutputComplexImpl(~)
    varargout = {true};
  end 
  
  function getInfoParameters(obj)
   
    % Setup parameters
    obj.pNt = obj.NumTransmitAntennas; 
    obj.pNr = obj.NumReceiveAntennas; 
    obj.pNumLinks = obj.pNt * obj.pNr;
    
    % Doppler components
    wavelength   = 3e8/obj.CarrierFrequency;                
    ScatterSpeed = obj.ScattererSpeed; % Scatterer speed in m/s;
    
    % Cut-off frequency f_D, in Hz
    obj.pDoppler = ScatterSpeed/wavelength;

    % Normalized Doppler spread f_D
    normalizationFactor = 1/300;

    % Channel sampling frequency
    obj.pFc = obj.pDoppler/normalizationFactor;
    
    targetOverSampleFactor = 10;
   
    % Limitation between sample rate and cut-off frequency factors
    coder.internal.errorIf(any(obj.SampleRate < ...
         targetOverSampleFactor * obj.pFc), ...
        'wlan:wlanChannel:MaxDopplerAndInputSamplingRate');
    
    coder.extrinsic('wlan.internal.spatialCorrelation'); 
    modelConfig = struct( ...
                  'NumTransmitAntennas',    obj.NumTransmitAntennas, ...
                  'NumReceiveAntennas',     obj.NumReceiveAntennas, ...
                  'TransmitAntennaSpacing', obj.TransmitAntennaSpacing, ...
                  'ReceiveAntennaSpacing',  obj.ReceiveAntennaSpacing, ...
                  'DelayProfile',           obj.DelayProfile, ...
                  'UserIndex',              obj.UserIndex, ...
                  'ChannelBandwidth',       obj.ChannelBandwidth, ...
                  'TransmitReceiveDistance',obj.TransmitReceiveDistance, ...
                  'CarrierFrequency',       obj.CarrierFrequency, ...
                  'TransmissionDirection',  obj.TransmissionDirection);
             
    if coder.target('MATLAB')
        out = wlan.internal.spatialCorrelation(modelConfig);                                 
    else
        out = coder.const((wlan.internal.spatialCorrelation(modelConfig)));
    end
    
    obj.pAntennaCorrelation  = out.AntennaCorrelation;
    obj.pPathPower           = out.PathPower;
    obj.pPathDelays          = out.PathDelays;
    obj.pPathloss            = out.Pathloss;
    obj.pDefaultShadowFading = out.ShadowFading;
    obj.pKfactor             = out.Kfactor;
    
    obj.pNumPaths = size(obj.pPathPower,2);
    
    % This is to display the initial shadow fading value through an info
    % function.
    obj.pShadowFading = out.ShadowFading;
    
    % Path power in dBs returned by the info function
    obj.pPathPowerdBs = 10*log10(obj.pPathPower);
    
    % Setup channel filter
    setupChannelFilter(obj);
      
  end
end

methods(Access = private)    
  
  function setupChannelFilter(obj)
  
    % Initial estimate of tapidx range. Minimum value of range is 0.
    err1 = 0.1;         % Small value
    c = 1/(pi*err1);    % Based on bound sinc(x) < 1/(pi*x)

    tRatio = (obj.pPathDelays * obj.SampleRate).';
    tapidx = min(floor(min(tRatio) - c), 0) : ceil(max(tRatio) + c);
    
    % Pre-compute channel filter tap gain sinc interpolation matrix
    if isempty(coder.target)   
        A1 = bsxfun(@minus, tRatio, tapidx);
        A = sinc(A1);
    else % Written in this way to avoid constant folding timeout error
        A1 = (repmat(tRatio, size(tapidx)));
        A  = (double(sinc(A1 - repmat(tapidx, size(tRatio)))));
    end
    
    % The following steps ensure that the tap index vector is shortened
    % when tRatio values are close to integer values.
    err2 = 0.01;
    if isempty(coder.target)   
        maxA = max(abs(A), [], 1); 
        significantIdx = find(maxA > (err2 * max(maxA)));
    else 
        maxA = (max(abs(A), [], 1)); 
        significantIdx = (double(find(maxA > (err2 * max(maxA)))));
    end
    
    t1 = min(tapidx(significantIdx(1)), 0);
    t2 = tapidx(significantIdx(end));

    if (t2 - t1) > 200
        coder.internal.warning('wlan:wlanChannel:PathDelayTooLarge');
    end

    % If the first tap index is negative, then the channel filter delay is
    % positive. This is the usual case. But if the first tap index is
    % positive, the channel filter delay is *negative*. This is an unusual
    % case for which the smallest path delay is much greater than the
    % input signal's sample period.
    obj.pChannelFilterDelay = round(tRatio(1)) - t1; 
    if obj.pChannelFilterDelay < 0
        coder.internal.warning('wlan:wlanChannel:ChannelFilterDelayNeg');
    end

    % Set up channel filter tap gain sinc interpolation matrix, [NP, nTaps]
    tapidx2 = t1:t2;
    if isempty(coder.target)
        A2 = bsxfun(@minus, tRatio, tapidx2);
        obj.pSincInterpMatrix = sinc(A2); 
    else % Written in this way to avoid constant folding timeout error
        A2 = (repmat(tRatio, size(tapidx2)));
        obj.pSincInterpMatrix = (double(sinc(A2 - repmat(tapidx2, size(tRatio)))));
    end
    
    obj.pNumChannelTaps = size(obj.pSincInterpMatrix, 2);  
  end
 
  function y = channelFilter(obj, x, fadingSamples)
    Ns = size(x, 1);    % Number of samples
    nTx = obj.pNt;      % Number of Tx
    nRx = obj.pNr;      % Number of Rx
    NP = obj.pNumPaths; % Number of paths
    
    % Initialize number of active links 
    activeTxIdx = 1:nTx;
    activeRxIdx = 1:nRx;
    activeLIdx  = 1:nTx*nRx;
    activeMIdx  = 1:nTx*nRx*NP;    
    
    numActiveTx = length(activeTxIdx);
    numActiveRx = length(activeRxIdx);
     
    if obj.pNumChannelTaps == 1 % Frequency-flat fading
        numActiveL  = length(activeLIdx);
        z = reshape(fadingSamples(:, activeMIdx), [Ns*numActiveL, NP]); 
        g = reshape(z * obj.pSincInterpMatrix, Ns, numActiveRx, []); 
        y = sum(bsxfun(@times, reshape(x, Ns, 1, []), g), 3);
        
    else % Frequency-selective fading
        filterState = obj.pChannelFilterState;
        %Calculate tap gain from path gain input
        % Note that conv((A*B), C) = A*conv(B, C). When B is a constant
        % matrix, right hand equation can have much faster
        % implementation than left-hand one by using 'filter' function.
         filterOut = coder.nullcopy(complex(zeros(Ns, numActiveTx*NP)));
         for i = 1:NP
            colIdx = (i-1)*numActiveTx + (1:numActiveTx);
            stateIdx = (i-1)*nTx + activeTxIdx;
            [filterOut(:,colIdx), filterState(:,stateIdx)] = ...
               filter(obj.pSincInterpMatrix(i, :), 1, ...
                          x(:,1:numActiveTx), filterState(:,stateIdx), 1);
         end
         y = sum(bsxfun(@times, reshape(fadingSamples(:,activeMIdx), ...
             Ns, numActiveRx, numActiveTx*NP), reshape(filterOut, ...
             Ns, 1, numActiveTx*NP)), 3);
         
         % Save filter state
         obj.pChannelFilterState = filterState;
    end 
    
  end
  
  function y = initializeRice(cfg)
                   
    nTx = cfg.NumTransmitAntennas;
    txSpacingNormalized = cfg.TransmitAntennaSpacing;
    txAoDLOSrad = pi/4;
    nRx = cfg.NumReceiveAntennas;
    rxSpacingNormalized = cfg.ReceiveAntennaSpacing;
    rxAoDLOSrad = pi/4;

    stepTx   = exp(1i*2*pi*txSpacingNormalized*sin(txAoDLOSrad));
    vectorTx = stepTx.^(0:nTx-1);
    stepRx   = exp(1i*2*pi*rxSpacingNormalized*sin(rxAoDLOSrad));
    vectorRx = stepRx.^((0:nRx-1).');
    y = vectorRx * vectorTx;
  end
  
  function [ricianLOS, ricianNLOS, pathPower] = ricianComponents(obj)
    
    % Computation of Rician steering matrix
    ricianMatrix = initializeRice(obj);
    
    % Initialization of the LOS component
    riceFactor = 10.^(.1.*obj.pKfactor);
    ricianLOS  = complex(zeros(obj.pNumPaths*obj.pNumLinks ,1));
    ricianNLOS = complex(zeros(obj.pNumPaths*obj.pNumLinks ,1));

    % Computation of the power delay profile of the (LOS+NLOS) power

    % The PDP is defined as the time dispersion of the NLOS power. The
    % addition of the LOS component modifies the time dispersion of the
    % total power.

    pathPower =obj.pPathPower.*(1+riceFactor);
        
    % Normalization of the power delay profile of the (LOS+NLOS) power
    if obj.NormalizePathGains
        pathPower = pathPower./sum(pathPower);
    end

    out = reshape(ricianMatrix, obj.pNumLinks, 1);

    for i = 1:obj.pNumPaths
        ricianLOS(1+(i-1)*obj.pNumLinks:i*obj.pNumLinks,1) = ...
                    sqrt(complex(pathPower(1,i))).* ...
                    sqrt(riceFactor(i)/(riceFactor(i)+1)).*out;
        ricianNLOS(1+(i-1)*obj.pNumLinks:i*obj.pNumLinks,1)= ...
                    sqrt(1/(riceFactor(i)+1)).*ones(obj.pNumLinks, 1);
    end;
  end
 
  function out = generateNoiseSamples(obj, N)
        
    NP = obj.pNumPaths;
    NL = obj.pNt * obj.pNr;
    M  = NP * NL;
    out = complex(zeros(NL*NP,N));
    const = sqrt(.5);
    
    % Gaussian noise generation
    if coder.target('MATLAB') 
        if strcmp(obj.RandomStream, 'Global stream')
            if obj.pLegacyGenerator
               stream = RandStream.getGlobalStream;
               stream.State = obj.pWGNState;
               wgnoise = const*complex(randn(stream, M, N), ...
                       randn(stream, M, N));
               obj.pWGNState = stream.State;
            else 
                % Separating legacy global syntax from this implementation
                w = (randn(2*M, N)).';
                wgnoise = const*(w(:,1:M) + 1i*w(:,M+1:end));
            end
        else 
            obj.pStream.State = obj.pWGNState; % Retrieve previous state
            if obj.pLegacyGenerator 
                wgnoise =  const*complex(randn(obj.pStream, M, N), ...
                           randn(obj.pStream, M, N));
            else
                w = (randn(obj.pStream, 2*M, N)).';
                wgnoise = const*(w(:,1:M) + 1i*w(:,M+1:end));
            end
            obj.pWGNState = obj.pStream.State; % Log randstream state
        end
    else
        if strcmp(obj.RandomStream, 'Global stream')
            if obj.pLegacyGenerator
               wgnoise = const*complex(randn(M, N), randn(M, N));
            else
               w =  (randn(2*M, N)).';
               wgnoise = const*(w(:,1:M) + 1i*w(:,M+1:end));
            end
        else
            wgnoise = coder.nullcopy(complex(zeros(M,N))); %#ok<NASGU>
            if obj.pLegacyGenerator 
		        w = coder.nullcopy(zeros(M,2*N));
                state = obj.pWGNState;
                for colIdx = 1:2*N 
                    for rowIdx= 1:M
                        [state, w(rowIdx,colIdx)] =  ...
                            eml_rand_mt19937ar('generate_normal', state);
                    end
                end
                wgnoise =  const*complex(w(:,1:N), w(:,N+1:end));
            else
                w = coder.nullcopy(zeros(N, 2*M));
                state = obj.pWGNState;
                for rowIdx = 1:N % Noise is generated column-wise
                    for colIdx = 1:2*M
                        [state, w(rowIdx, colIdx)] = eml_rand_mt19937ar('generate_normal', state);
                    end
                end
                 wgnoise = const*(w(:,1:M) + 1i*w(:,M+1:end)); 
            end
            obj.pWGNState = state;
        end
    end
    
    % Filter definition
    filterOrder = 7;

    switch obj.DelayProfile

         case 'Model-F'
            % Indoor, bell shape Doppler spectrum with spike
            filterStateIn = obj.pFilterState;

            obj.pFilterState = complex(zeros(filterOrder, ...
                                    obj.pNt*obj.pNr*obj.pNumPaths));
            SpikePos = 3; %tap with Bell+Spike shape
            %index of Bell shaped fading coefficients
            BellTapsPositions = [1:(SpikePos-1)*NL ((SpikePos)*NL+1):M]; 
            %index of Bell+Spike shaped fading coefficients
            SpikeTapsPositions = ((SpikePos-1)*NL+1):(SpikePos)*NL;
            
            if obj.pLegacyGenerator
                %Bell shaped fading
                [FadingMatrixTimeBellTap,  ...
                    obj.pFilterState(:,BellTapsPositions)] =...
                    filter(obj.pNumeratorB, obj.pDenominatorB,  ...
                    wgnoise(BellTapsPositions,:), ...
                    filterStateIn(:,BellTapsPositions), 2);

                %Bell+Spike shaped fading
                [FadingMatrixTimeSpikeTap, ...
                    obj.pFilterState(:,SpikeTapsPositions)] = ...
                    filter(obj.pNumeratorS, obj.pDenominatorS, ...
                    wgnoise(SpikeTapsPositions,:), ...
                    filterStateIn(:, SpikeTapsPositions), 2);
                            
                % Fading Matrix Time with bell+spike on the 3-rd tap
                out(BellTapsPositions,:) = FadingMatrixTimeBellTap;
                out(SpikeTapsPositions,:) = FadingMatrixTimeSpikeTap;               
            else
                % Bell shaped fading
                [FadingMatrixTimeBellTap,  ...
                    obj.pFilterState(:,BellTapsPositions)] =...
                    filter(obj.pNumeratorB, obj.pDenominatorB,  ...
                    wgnoise(:,BellTapsPositions), ...
                    filterStateIn(:,BellTapsPositions), 1);

                %Bell+Spike shaped fading
                [FadingMatrixTimeSpikeTap, ...
                    obj.pFilterState(:,SpikeTapsPositions)] = ...
                    filter(obj.pNumeratorS, obj.pDenominatorS, ...
                    wgnoise(:,SpikeTapsPositions), ...
                    filterStateIn(:,SpikeTapsPositions), 1);
                
                % Fading Matrix Time with bell+spike on the 3-rd tap
                out(BellTapsPositions,:) = FadingMatrixTimeBellTap.';
                out(SpikeTapsPositions,:) = FadingMatrixTimeSpikeTap.';    
            end
        otherwise   
            % Indoor, bell shape Doppler spectrum
            if obj.pLegacyGenerator
                [out,obj.pFilterState] = filter(obj.pNumeratorB,obj.pDenominatorB, ...
                                  wgnoise, obj.pFilterState, 2);

            else
                [filterOut,obj.pFilterState] = filter(obj.pNumeratorB, obj.pDenominatorB, ...
                                            wgnoise, obj.pFilterState, 1);
                out =  filterOut.';  
                  
            end
    end

  end

 % Add fluorescent effects
 function out =  fluorescentEffects(obj, FadingMatrix, fadingSamplingTime)
        
    % Amplitudes of the modulating signal
    amplitudes = [0 -15 -20];
    % Taps to be used in model D
    tapsModelD = [12 14 16];
    % Taps to be used in model E
    tapsModelE = [3 5 7];

    % Initialization of the modulation matrix
    ModulationMatrix = complex(zeros(size(FadingMatrix)));

    % Computation of the modulation function g(t)
    g = complex(zeros(3,size(fadingSamplingTime,2)));

    for k = 1:3
        g(k,:) = (10^(amplitudes(k)/20)).* ...
            exp(1i.*(4.*pi.*(2*(k-1)+1).* ...
            obj.pPowerLineFrequency.*fadingSamplingTime ...
            + obj.pRandomPhase));
    end;

    modulationFunction = sum(g);

    % Modulation matrix
    if strcmpi(obj.DelayProfile, 'Model-D')
        
        for ii=1:obj.pNt
            for jj=1:obj.pNr
                for k=1:3
                    ModulationMatrix((((tapsModelD(k)-1)*obj.pNt*obj.pNr) ...
                        +((ii-1)*obj.pNr)+jj),:) = modulationFunction;
                end;
            end;   
        end;

    else       
        for ii=1:obj.pNt
            for jj=1:obj.pNr
                for k=1:3
                    ModulationMatrix((((tapsModelE(k)-1)*obj.pNt*obj.pNr) ...
                        +((ii-1)*obj.pNr)+jj),:) = modulationFunction;
                end;
            end;   
        end;
    end;

    % Multiplication of the modulation matrix and the fading matrix
    ModulationMatrix = ModulationMatrix.*FadingMatrix;
    
    if obj.pLegacyGenerator
        % Calculation of the energy of both matrices
        ModulationMatrixEnergy = sum(sum(abs(ModulationMatrix).^2));
        FadingMatrixEnergy     = sum(sum(abs(FadingMatrix).^2));

        % Comparison and calculation of the normalization constant alpha
        RealInterfererToCarrierRatio = ModulationMatrixEnergy/FadingMatrixEnergy;
        NormalizationConstant = sqrt(double( ...
            obj.pRandomInterfererToCarrierRatio/RealInterfererToCarrierRatio));

        out = bsxfun(@plus, FadingMatrix, NormalizationConstant(1,1).* ModulationMatrix);
    else
        % Calculation of the energy of both matrices
        ModulationMatrixEnergy = (sum(abs(ModulationMatrix).^2));
        FadingMatrixEnergy     = (sum(abs(FadingMatrix).^2));

        % Comparison and calculation of the normalization constant alpha
        RealInterfererToCarrierRatio = ModulationMatrixEnergy./FadingMatrixEnergy;
        NormalizationConstant = sqrt(double( ...
            obj.pRandomInterfererToCarrierRatio./RealInterfererToCarrierRatio));

        NC = repmat(NormalizationConstant, size(ModulationMatrix,1),1);

        out = bsxfun(@plus, FadingMatrix, NC.* ModulationMatrix);
    end

 end
 
 function out = generateWGNoise(obj)
   if coder.target('MATLAB')   
      if strcmp(obj.RandomStream, 'Global stream')
           if obj.pLegacyGenerator 
               obj.pWGNState = obj.pStream.State;
               out = randn(obj.pStream,1,1); 
           else
               out = randn(1,1);
           end
      else
        % Continuous generation of samples unlike TGn
        obj.pStream.State = obj.pWGNState;
        out =  randn(obj.pStream,1,1);
        obj.pWGNState = obj.pStream.State;
      end
    else
      if strcmp(obj.RandomStream, 'Global stream')
         out = randn(1,1);
      else 
         [obj.pWGNState, out] = eml_rand_mt19937ar('generate_normal', obj.pWGNState);
     end
   end
 end
 
 function H = generateFadingSamples(obj,runningTime,numChSamples)
     
       fadingSamplingTime = 1/obj.pFc;
       nTx = obj.pNt;
       nRx = obj.pNr;
     
       % Time references
       obj.pFadingTime  = obj.pFadingOffset-(obj.pFadingOffset~= 0)*2*fadingSamplingTime: fadingSamplingTime:...
                     obj.pFadingOffset+((numChSamples-1)*fadingSamplingTime); 

       obj.pOversampledFadingTime = obj.pFadingTime(1):.1*fadingSamplingTime:...
                                  obj.pFadingTime(size(obj.pFadingTime,2));
        
                              
       % Computation of the matrix of fading coefficients 
        newFadingMatrixTime = generateNoiseSamples(obj, numChSamples);
        
       % State retention
       if obj.pOffset==0;
           fadingMatrixTime = newFadingMatrixTime;
       else
           fadingMatrixTime = [obj.pLastTwoColumn, newFadingMatrixTime]; 
       end
        
       obj.pLastTwoColumn = newFadingMatrixTime(:, ...
                                        numChSamples-1:numChSamples);

       % Spatial correlation of the fading matrix
       fadingMatrixTime = obj.pAntennaCorrelation * fadingMatrixTime;

       pdpCoef =kron(sqrt(complex(obj.pPathPower)), ones(1, nTx*nRx));
 
       % Normalization of the correlated fading processes
       fadingMatrixTime = diag(pdpCoef)*fadingMatrixTime;

       % Calculation of the Rice phasor AoA/AoD hard-coded to 45 degrees
       ricePhasor = exp(1i.*2.*pi.*obj.pDoppler.*cos(pi/4).*obj.pFadingTime);

       % Addition of the Rice component
       fadingMatrix = (obj.pRicianLOS*ricePhasor)+(obj.pRicianNLOS ...
                      *ones(1,size(fadingMatrixTime,2))).*fadingMatrixTime;
 
       % First (partial) interpolation, at 10*FadingSamplingFrequency_Hz
       % to avoid aliasing of the fluorescent effect
       obj.pOversampledFadingMatrix = (interp1(obj.pFadingTime.',fadingMatrix.', ...
                                obj.pOversampledFadingTime.','linear')).';
 
       % Call to generate fluorescent light effects in models D and E
       if ( strcmpi(obj.DelayProfile, 'Model-D') || ...
                                strcmpi(obj.DelayProfile, 'Model-E')) && obj.FluorescentEffect
            obj.pOversampledFadingMatrix = fluorescentEffects(obj,  ...
                obj.pOversampledFadingMatrix, obj.pOversampledFadingTime);
       end;

       % Large-scale fading
       if ~strcmp(obj.LargeScaleFadingEffect,'None')
            obj.pOversampledFadingMatrix =  ...
                obj.pOversampledFadingMatrix.*obj.pLargeScaleLosses;
       end

       % Second interpolation, at SamplingRate_Hz
       H = interp1(obj.pOversampledFadingTime.',obj.pOversampledFadingMatrix.', ...
                  runningTime.','linear');
     
 end
 
end

end
 
% [EOF]
