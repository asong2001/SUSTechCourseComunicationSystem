function [bits, eqDataSym, varargout] = wlanVHTSIGBRecover(rxVHTSIGB, ...
    chanEst, noiseVarEst, chanBW, varargin)
%WLANVHTSIGBRECOVER Recover information bits in VHT-SIG-B field
% 
%   BITS = wlanVHTSIGBRecover(RXVHTSIGB, CHANEST, NOISEVAREST, CHANBW)
%   recovers the information bits in the VHT-SIG-B field.
%
%   BITS is an int8 column vector containing the recovered information
%   bits.
%
%   RXVHTSIGB is the received time-domain VHT-SIG-B field signal. It is a
%   Ns-by-Nr matrix of real or complex values, where Ns represents the
%   number of time-domain samples in the VHT-SIG-B field and Nr represents
%   the number of receive antennas. Ns can be greater than the VHT-SIG-B
%   field length; in this case additional samples at the end of RXVHTSIGB
%   are not used.
% 
%   CHANEST is the estimated channel at data and pilot subcarriers based on
%   the VHT-LTF. It is a real or complex array of size Nst-by-Nsts-by-Nr,
%   where Nst represents the total number of occupied subcarriers.
% 
%   NOISEVAREST is the noise variance estimate. It is a real, nonnegative
%   scalar.
%  
%   CHANBW is the channel bandwidth, one of 'CBW20', 'CBW40', 'CBW80' and
%   'CBW160'.
% 
%   BITS = wlanVHTSIGBRecover(..., CFGREC) allows different algorithm
%   options for information bit recovery via the input CFGREC, which is a
%   <a href="matlab:help('wlanRecoveryConfig')">wlanRecoveryConfig</a> configuration object. When the CFGREC input is not  
%   specified, the default property values of the <a href="matlab:help('wlanRecoveryConfig')">wlanRecoveryConfig</a>  
%   object are adopted in the recovery.
%
%   [..., EQDATASYM, CPE] = wlanVHTSIGBRecover(...) also returns the
%   equalized subcarriers and common phase error.
%
%   EQDATASYM is a complex column vector of length Nsd containing the
%   equalized symbols at data subcarriers. Nsd represents the number of
%   data subcarriers.
%
%   CPE is a scalar containing the common phase error between the received
%   and expected OFDM symbol.
% 
%   Example: 
%   %  Recover the information bits in VHT-SIG-B field via channel 
%   %  estimation on VHT-LTF over a 4 x 2 quasi-static fading channel
%
%     % Configure a VHT configuration object 
%     chanBW = 'CBW80';
%     cfgVHT = wlanVHTConfig( 'ChannelBandwidth',    chanBW, ...
%         'NumTransmitAntennas', 4, 'NumSpaceTimeStreams', 4); 
%  
%     % Generate VHT-LTF and VHT-SIG-B field signals
%     txVHTLTF = wlanVHTLTF(cfgVHT); 
%     [txVHTSIGB, VHTSIGBBits] = wlanVHTSIGB(cfgVHT);
% 
%     % Pass through a 4 x 2 quasi-static fading channel with AWGN 
%     H = 1/sqrt(2)*complex(randn(4, 2), randn(4, 2));
%     rxVHTLTF  = awgn(txVHTLTF  * H, 10);
%     rxVHTSIGB = awgn(txVHTSIGB * H, 10);
% 
%     % Perform channel estimation based on VHT-LTF
%     demodVHTLTF = wlanVHTLTFDemodulate(rxVHTLTF, cfgVHT, 1);
%     chanEst = wlanVHTLTFChannelEstimate(demodVHTLTF, cfgVHT);
% 
%     % Recover information bits in VHT-SIG-B
%     recVHTSIGBBits = wlanVHTSIGBRecover(rxVHTSIGB, chanEst, 0.1, chanBW);
%
%     % Compare against original information bits
%     disp(isequal(VHTSIGBBits, recVHTSIGBBits));
% 
%   See also wlanRecoveryConfig, wlanVHTSIGB, wlanVHTLTF,
%   wlanVHTLTFDemodulate, wlanVHTLTFChannelEstimate.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

narginchk(4, 5);
nargoutchk(0, 3);

% Calculate CPE if requested
if nargout>2
    calculateCPE = true;
else
    calculateCPE = false;
end

% Validate channel bandwidth input
coder.internal.errorIf(~ischar(chanBW) || ...
    ~any(strcmp(chanBW, {'CBW20','CBW40','CBW80','CBW160'})), ...
    'wlan:wlanVHTSIGBRecover:InvalidChanBW');

numSTS = size(chanEst, 2);

% Get OFDM configuration
[cfgOFDM,dataInd,pilotInd] = wlanGetOFDMConfig(chanBW, 'Long', 'VHT', ...
    numSTS);

% Validate VHT-SIG-B field signal input
validateattributes(rxVHTSIGB, {'double'}, {'2d','finite','nrows', ...
    cfgOFDM.FFTLength*5/4}, 'rxVHTSIGB', 'VHT-SIG-B field signal'); 
numRx = size(rxVHTSIGB, 2);

% Validate channel estimates
validateattributes(chanEst, {'double'}, {'3d','finite','nonempty'}, ...
    'chanEst', 'channel estimates'); 

% Cross validate inputs
numST = numel([dataInd; pilotInd]); % Total number of occupied subcarriers
coder.internal.errorIf(size(chanEst, 1) ~= numST || ...
    (size(chanEst, 2) > 8) || (size(chanEst, 3) ~= numRx), ...
    'wlan:wlanVHTSIGBRecover:InvalidChanEst', numST, numRx);

% Extract data and pilot subcarriers from channel estimate
chanEstData = chanEst(dataInd,:,:);
chanEstPilots = chanEst(pilotInd,:,:);

% Validate noise variance estimate input
validateattributes(noiseVarEst, {'double'}, ...
    {'real','scalar','nonnegative','finite'}, ...
    'noiseVarEst', 'noise variance estimate'); 

% Validate optional recovery configuration input
if nargin == 5
    validateattributes(varargin{1}, {'wlanRecoveryConfig'}, {'scalar'}, ...
        'cfgRec', 'recovery configuration object');
    symOffset = varargin{1}.OFDMSymbolOffset;
    eqMethod  = varargin{1}.EqualizationMethod; 
    pilotPhaseTracking = varargin{1}.PilotPhaseTracking;
else
    symOffset = 0.75;
    eqMethod  = 'MMSE'; 
    pilotPhaseTracking = 'PreEQ';
end

% OFDM demodulation
[ofdmDemodData, ofdmDemodPilots] = ... % data is of size [numSD, 1, numRx]
    wlanOFDMDemodulate(rxVHTSIGB, cfgOFDM, symOffset); 

% Pilot phase tracking
if calculateCPE==true || strcmp(pilotPhaseTracking, 'PreEQ')
    % Get reference pilots, from Eqn 22-47, IEEE Std 802.11ac-2013
    z = 3; % Offset by 3 to allow for L-SIG and VHT-SIG-A pilot symbols
    refPilots = vhtPilots(1, z, chanBW, numSTS);
    
    % Estimate CPE and phase correct symbols
    cpe = commonPhaseErrorEstimate(ofdmDemodPilots, chanEstPilots, ...
        refPilots);
    if strcmp(pilotPhaseTracking, 'PreEQ')
        ofdmDemodData = commonPhaseErrorCorrect(ofdmDemodData, cpe);
    end
    varargout{1} = cpe.'; % Permute to Nsym-by-1
end

% Perform equalization
% Flip the 4th and 8th STS, i.e., P matrix multiplication
if any(numSTS == [4 7 8])
    chanEstData(:,4:4:end,:) = -chanEstData(:,4:4:end,:);
end
[eqDataSym, csiData] = wlanEqualize(ofdmDemodData, sum(chanEstData,2), ...
    eqMethod, noiseVarEst); % [numSD, 1]

% Constellation demapping: For VHT-SIG-B, avoid calling segment deparser
% which just converts a column vector into 2 columns.
qamDemodOut = wlanConstellationDemodulate(eqDataSym, 1, noiseVarEst);
qamDemodOut = qamDemodOut .* csiData; % [numSD, 1]

% Deinterleaving & segment parser for CBW160
num20 = cfgOFDM.FFTLength/64;
numSD = length(cfgOFDM.DataIndices);
if strcmp(chanBW, 'CBW160')
    deparserOut = zeros(size(qamDemodOut));
    deparserOut(1:2:end, :) = wlanBCCDeinterleave( ...
        qamDemodOut(1:end/2,:), 'VHT', numSD/2, 1, 20*num20, 1);
    deparserOut(2:2:end, :) = wlanBCCDeinterleave( ... 
        qamDemodOut(end/2+1:end,:), 'VHT', numSD/2, 1, 20*num20, 1);
else
    deparserOut = wlanBCCDeinterleave(qamDemodOut, 'VHT', numSD, 1, ...
                                      20*num20, 1);
end

% Remove redundant zeros between information bit repetitions
if strcmp(chanBW, 'CBW80')
    infoBitRep = deparserOut(1:end-2, :);
elseif strcmp(chanBW, 'CBW160')
    infoBitRep = deparserOut([1:end/2-2,end/2+1:end-2], :);
else
    infoBitRep = deparserOut;
end

% BCC decoding: length 26 for 'CBW20', 27 for 'CBW40' and 29 for 'CBW80'
% and 'CBW160'
recBitLen = length(infoBitRep)/num20/2;
bits = wlanBCCDecode(mean(reshape(infoBitRep, [], num20), 2), ...
        '1/2', recBitLen);

end

% [EOF]