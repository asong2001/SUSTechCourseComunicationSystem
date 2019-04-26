function est = wlanVHTLTFChannelEstimate(rxSym,cfgVHT,varargin)
% wlanVHTLTFChannelEstimate Channel estimation using the VHT-LTF
%   EST = wlanVHTLTFChannelEstimate(RXSYM,CFGVHT) returns the estimated
%   channel between all space-time streams and receive antennas using the
%   Very High Throughput Long Training Field (VHT-LTF). The channel
%   estimate includes the effect of the applied spatial mapping matrix and
%   cyclic shifts at the transmitter.
%
%   EST is a complex Nst-by-Nsts-by-Nr array containing the estimated
%   channel at data and pilot subcarriers, where Nst is the number of
%   occupied subcarriers, Nsts is the number of space-time streams and Nr
%   is the number of receive antennas.
%
%   RXSYM is a complex Nst-by-Nsym-by-Nr array containing demodulated
%   VHT-LTF OFDM symbols. Nsym is the number of demodulated VHT-LTF
%   symbols.
%
%   CFGVHT is a packet format configuration object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a>. 
%
%   EST = wlanVHTLTFChannelEstimate(...,SPAN) performs frequency smoothing
%   by using a moving average filter across adjacent subcarriers to reduce
%   the noise on the channel estimate. The span of the filter in
%   subcarriers, SPAN, must be odd. If adjacent subcarriers are highly
%   correlated frequency smoothing will result in significant noise
%   reduction, however in a highly frequency selective channel smoothing
%   may degrade the quality of the channel estimate.
%
%   Examples:
%
%    Example 1:
%    %  Generate a time domain waveform txWaveform for an 802.11ac VHT
%    %  packet and combine all transmit antennas onto one receive antenna.
%    %  Extract and demodulate the VHT-LTF and perform channel estimation
%    %  with a frequency smoothing span of 3.
%
%       cfgVHT = wlanVHTConfig;         % Create packet configuration
%       cfgVHT.NumTransmitAntennas = 2; % 2 transmit antennas
%       cfgVHT.NumSpaceTimeStreams = 2; % 2 spatial streams
%       txWaveform = wlanWaveformGenerator([1;0;0;1],cfgVHT);
%
%       rxWaveform = sum(txWaveform,2); % Combine all transmit antennas
%       idnVHTLTF = wlanFieldIndices(cfgVHT,'VHT-LTF');
%       sym = wlanVHTLTFDemodulate(rxWaveform(idnVHTLTF(1):idnVHTLTF(2),:), ...
%                                   cfgVHT);
%       smoothingSpan = 3;
%       est = wlanVHTLTFChannelEstimate(sym,cfgVHT,smoothingSpan);
%
%    Example 2:
%    %  Generate a time domain waveform for an 802.11ac VHT packet, pass it
%    %  through a TGac fading channel and perform VHT-LTF channel 
%    %  estimation.
%
%       cfgVHT = wlanVHTConfig;       % Create packet configuration
%       txWaveform = wlanWaveformGenerator([1;0;0;1],cfgVHT);
% 
%       % Configure channel
%       tgac = wlanTGacChannel;
%       tgac.SampleRate = 80e6;
%       tgac.ChannelBandwidth = 'CBW80';
% 
%       % Pass through channel (with zeros to allow for channel delay)
%       rxWaveform = step(tgac,[txWaveform; zeros(15,1)]);
%       rxWaveform = rxWaveform(5:end,:); % Synchronize for channel delay
% 
%       % Extract VHT-LTF and perform channel estimation
%       indVHTLTF = wlanFieldIndices(cfgVHT,'VHT-LTF');
%       sym = wlanVHTLTFDemodulate(rxWaveform( ...
%                indVHTLTF(1):indVHTLTF(2),:),cfgVHT);
%       est = wlanVHTLTFChannelEstimate(sym,cfgVHT);
%
%   See also wlanVHTConfig, wlanVHTLTFDemodulate, wlanVHTDataRecover,
%            wlanLLTFChannelEstimate, wlanHTLTFChannelEstimate.
 
%   Copyright 2015 The MathWorks, Inc.

%#codegen

% Validate number of arguments
narginchk(2,3);

if nargin > 2
    span = varargin{1};
    enableFreqSmoothing = true;
else
    % Default no frequency smoothing
    enableFreqSmoothing = false;
end

% Validate the packet format configuration object is a valid type
validateattributes(cfgVHT,{'wlanVHTConfig'},{'scalar'}, ...
    'wlanVHTLTFChannelEstimate','packet format configuration object');
coder.internal.errorIf(cfgVHT.NumUsers ~= 1, ...
    'wlan:wlanChannelEstimate:SingleUserOnly');

% Validate symbol type
validateattributes(rxSym,{'single','double'},{'3d'}, ...
    'wlanVHTLTFChannelEstimate','VHT-LTF OFDM symbol(s)');

cbw = cfgVHT.ChannelBandwidth;
numSC = size(rxSym,1);
numRxAnts = size(rxSym,3);
numSTS = cfgVHT.NumSpaceTimeStreams;

% Return an empty if empty symbols
if isempty(rxSym)
    est = zeros(numSC,numSTS,numRxAnts);
    return;
end

% Verify number of subcarriers to estimate
[cfgOFDM,dataInd] = wlanGetOFDMConfig(cbw,'Long','VHT',numSTS); % Get OFDM configuration
FFTLen = cfgOFDM.FFTLength;
[ind,sortIdx] = sort([cfgOFDM.DataIndices; cfgOFDM.PilotIndices]);
coder.internal.errorIf(numSC~=numel(ind), ...
    'wlan:wlanChannelEstimate:IncorrectNumSC',numel(ind),numSC);

% Get cyclic shifts applied at transmitter
csh = getCyclicShiftVal('VHT',numSTS,FFTLen/64*20);

if numSTS==1
    % Perform channel estimation for all subcarriers
    est = vhtltfEstimate(rxSym,cbw,numSTS,ind); 
else
    % Perform channel estimation for data carrying subcarriers as we
    % must interpolate the pilots
    estData = vhtltfEstimate(rxSym(dataInd,:,:),cbw,numSTS,cfgOFDM.DataIndices);

    % Undo cyclic shift for each STS before averaging and interpolation
    estData = wlanCyclicShift(estData,csh,FFTLen,'Rx',cfgOFDM.DataIndices);

    % Estimate pilot subcarriers
    estPilots = pilotInterpolation(estData,FFTLen,cfgOFDM.DataIndices, ...
        cfgOFDM.PilotIndices);

    % Combine data and pilots into one container
    allSubcarriers = [estData; estPilots];
    est = allSubcarriers(sortIdx,:,:);              
end

% Perform frequency smoothing
if enableFreqSmoothing
    % Smooth segments between DC gaps
    switch cbw
        case 'CBW20'
            numGroups = 1;
        case 'CBW40'
            numGroups = 2;
        case 'CBW80'
            numGroups = 2;
        otherwise % 'CBW160'
            numGroups = 4;
    end
    groupSize = size(est,1)/numGroups;
    for i = 1:numGroups
        idx = (1:groupSize)+(i-1)*groupSize;
        est(idx,:,:) = frequencySmoothing(est(idx,:,:),span);
    end
end  

% Re-apply cyclic shift after averaging and interpolation
if numSTS>1  
    est = wlanCyclicShift(est,csh,FFTLen,'Tx',ind);
end
end

%--------------------------------------------------------------------------
function estPilots = pilotInterpolation(estData,Nfft,dataIndices,pilotIndices)
% Interpolate over the pilot locations
    numSTS = size(estData,2);
    numRxAnts = size(estData,3);

    % Construct full FFT size to allow us to interpolate over DC nulls
    est = complex(ones(Nfft,numSTS,numRxAnts),ones(Nfft,numSTS,numRxAnts));
    est(dataIndices,:,:) = estData;
       
    % Interpolate over missing parts of the waveform in magnitude and
    % phase (as opposed to real and imaginary)
    magPart = interp1(dataIndices,abs(est(dataIndices,:,:)),1:Nfft);
    phasePart = interp1(dataIndices,unwrap(angle(est(dataIndices,:,:))),1:Nfft);
    [realPart,imagPart] = pol2cart(phasePart,magPart);
    estInterp = complex(realPart,imagPart);
    if isrow(estInterp)
        est = estInterp(:,:,1).';
    else
        est = estInterp;
    end
    
    % Extract pilots
    estPilots = est(pilotIndices,:,:);
end