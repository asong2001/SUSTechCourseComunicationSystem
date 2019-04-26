function est = helperNoiseEstimate(rxSym,varargin)
%helperNoiseEstimate Estimate noise power using L-LTF 
%
%   EST = HELPERNOISEESTIMATE(rxSym) estimates the mean noise power in
%   watts using the demodulated L-LTF symbols, assuming 1ohm resistance.
%   The estimated noise power in Non-HT packet is averaged over the number
%   of receive antennas.
%
%   RXSYM is the frequency-domain signal corresponding to the L-LTF.
%   It is a complex matrix or 3-D array of size Nst-by-2-by-Nr, where Nst
%   represents the number of used subcarriers in the L-LTF, and Nr
%   represents the number of receive antennas. Two OFDM symbols in the
%   L-LTF field are used to estimate the noise power.
%
%   EST = HELPERNOISEESTIMATE(RXSYM,CFGFORMAT) returns the estimated noise
%   power in VHT or HT fields using CFGFORMAT. The CFGFORMAT is a
%   format configuration object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a>, or <a
%   href="matlab:help('wlanHTConfig')">wlanHTConfig</a>. 
%   The number of subcarriers used within each field and the scaling
%   applied during demodulation differs between the Non-HT and HT/VHT
%   fields. Therefore the estimated noise power after demodulation in
%   HT/VHT fields is calculated by scaling the estimated noise power in the
%   L-LTF.
%
%   EST = HELPERNOISEESTIMATE(...,'Per Antenna') specifies the option of
%   estimating the noise for each receive antenna. When this option is
%   specified EST is a row vector of length Nr.
%
%   Example:
%   %  Estimate the noise variance of an HT packet.
%
%      cfgHT = wlanHTConfig;   % Create packet configuration
%      noisePower = -20;
%      AWGN = comm.AWGNChannel;
%      AWGN.NoiseMethod = 'Variance';
%      AWGN.Variance = 10^(noisePower/10);
%
%      Nst = 56;  % Data and pilot OFDM subcarriers in 20MHz, HT format
%      Nfft = 64; % FFT size for 20MHz bandwidth
%      nVarHT = 10^(noisePower/10)*(Nst/Nfft); % Non-HT noise variance
%
%      NumRxAnts = 1;
%      % Average noise estimate over 100 independent noise realization
%      for n=1:100
%         % Generate LLTF and add noise
%         rxSym = step(AWGN,wlanLLTF(cfgHT));
%         y = wlanLLTFDemodulate(rxSym,cfgHT);
%         noiseEst(n) = helperNoiseEstimate(y,cfgHT);
%      end
%       
%      % Check noise variance estimates without Channel 
%      noiseEstError = 10*log10(mean(noiseEst))-10*log10(nVarHT);
%      disp(['Error between noise variance and mean estimated noise ', ...
%      'power(dB): ' num2str(noiseEstError,'%2.2f ')]);
%
%   See also wlanLLTF, wlanLLTFDemodulate.
 
%   Copyright 2015 The MathWorks, Inc.

%#codegen

narginchk(1,3);

% Validate symbol type
validateattributes(rxSym,{'double'},{'3d','finite'}, ...
    'helperNoiseEstimate','L-LTF OFDM symbol(s)');

% Two L-LTF symbols are required to estimate the noise
coder.internal.errorIf(size(rxSym,2)~=2, ...
        'wlan:helperNoiseEstimate:IncorrectNumSyms');

numSC = size(rxSym,1);

if nargin == 2
    % helperNoiseEstimate(rxSym,'Per Antenna')
    varType = varargin{1};
    if isa(varType,'char') % Second input is set to 'Per Antenna'
        coder.internal.errorIf( ~( strcmpi(varType,'Per Antenna')), ...
        'wlan:helperNoiseEstimate:InvalidNoiseEstType');
    
        average = false;
        scalingFactor = 1; % Noise scaling factor
    else 
        % helperNoiseEstimate(rxSym,cfgVHT)
        validateInput(varType,numSC); % Validate input
       
        average = true;
        scalingFactor = noiseScaling(varType,numSC); % Noise scaling factor
    end
elseif nargin == 3
    % helperNoiseEstimate(rxSym,cfgVHT,'Per Antenna')
    varType1 = varargin{1};

    validateInput(varType1,numSC); % Validate input
    
    varType2 = varargin{2};
    coder.internal.errorIf( ~( strcmpi(varType2,'Per Antenna')), ...
        'wlan:helperNoiseEstimate:InvalidNoiseEstType');
    
    average = false;
    % Noise scaling factor
    scalingFactor = noiseScaling(varType1,numSC);
else 
    % helperNoiseEstimate(rxSym)
    average = true;
    scalingFactor = 1;
end

% Noise estimate
noiseEst = sum(abs(rxSym(:,1,:) - rxSym(:,2,:)).^2,1)/(2*numSC);
if average
    noise = mean(noiseEst);
else
    noise = squeeze(noiseEst).';
end

% Scale
est = noise*scalingFactor; 

end

function out = noiseScaling(cfgFormat,numSC)
    % Get the number of occupied subcarriers in HT and VHT fields. The
    % number of used subcarriers for HT and NonHT are same therefore
    % fix the string input of the following helper function to VHT
    [vhtData,vhtPilots] = helperSubcarrierIndices(cfgFormat,'VHT');
    NstVHT = numel(vhtData)+numel(vhtPilots);
    out = (NstVHT/numSC)*cfgFormat.NumSpaceTimeStreams;
end

function validateInput(cfgFormat,numSC)

    % Validate format type to be HT and VHT
    validateattributes(cfgFormat, ...
    {'wlanVHTConfig','wlanHTConfig'},{},mfilename,'second argument');
    
    chBW = cfgFormat.ChannelBandwidth;
    % Get number of used subcarriers in NonHT format
    [nonhtData,nonhtPilots] = helperSubcarrierIndices(chBW,'Legacy');
    nonhtNst = numel(nonhtData)+numel(nonhtPilots);
    
    % Verify number of subcarriers
    coder.internal.errorIf(numSC~=nonhtNst, ...
        'wlan:helperNoiseEstimate:IncorrectNumSC',nonhtNst,numSC);
end

% EOF