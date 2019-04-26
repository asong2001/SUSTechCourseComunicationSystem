function [bits, failCheck, eqDataSym, varargout] = wlanLSIGRecover( ...
    rxLSIG, chanEst, noiseVarEst, chanBW, varargin)
%WLANLSIGRECOVER Recover information bits in L-SIG field
% 
%   [BITS, FAILCHECK] = wlanLSIGRecover(RXLSIG, CHANEST, NOISEVAREST,
%   CHANBW) recovers the information bits in the L-SIG field and performs
%   parity and rate checks.
% 
%   BITS is an int8 column vector of length 24 containing the recovered
%   information bits.
% 
%   FAILCHECK is a logical scalar which is true if BITS fails the parity
%   check or its first 4 bits is not one of the eight legitimate rates.
%
%   RXLSIG is the received time-domain L-SIG field signal. It is a Ns-by-Nr
%   matrix of real or complex values, where Ns represents the number of
%   time-domain samples in the L-SIG field and Nr represents the number of
%   receive antennas. Ns can be greater than the L-SIG field length; in
%   this case additional samples at the end of RXLSIG are not used.
% 
%   CHANEST is the estimated channel at data and pilot subcarriers based on
%   the L-LTF. It is a real or complex array of size Nst-by-1-by-Nr, where
%   Nst represents the total number of occupied subcarriers.  The singleton
%   dimension corresponds to the single transmitted stream in the L-LTF
%   which includes the combined cyclic shifts if multiple transmit antennas
%   are used.
%   
%   NOISEVAREST is the noise variance estimate. It is a real nonnegative
%   scalar.
% 
%   CHANBW is a string describing the channel bandwidth and must one of:
%   'CBW5','CBW10','CBW20','CBW40','CBW80','CBW160'.
% 
%   [BITS, FAILCHECK] = wlanLSIGRecover(..., CFGREC) allows different
%   algorithm options for information bit recovery via the input CFGREC,
%   which is a <a href="matlab:help('wlanRecoveryConfig')">wlanRecoveryConfig</a> configuration object. When the CFGREC 
%   input is not specified, the default property values of the
%   <a href="matlab:help('wlanRecoveryConfig')">wlanRecoveryConfig</a> object are adopted in the recovery.
%
%   [..., EQDATASYM, CPE] = wlanLSIGRecover(...) also returns the equalized
%   subcarriers and common phase error.
%
%   EQDATASYM is a complex column vector of length 48 containing the
%   equalized symbols at data carrying subcarriers. There are 48 data
%   carrying subcarriers in the L-SIG field.
%
%   CPE is an scalar containing the common phase error between the received
%   and expected OFDM symbol.
%
%   Example: 
%   % Recover the information bits in L-SIG via channel estimation on L-LTF
%   % over a 3 x 2 quasi-static fading channel.
%
%     % Configure a VHT configuration object 
%     chanBW = 'CBW40'; 
%     cfgVHT = wlanVHTConfig( ...
%         'ChannelBandwidth', chanBW, ...
%         'NumTransmitAntennas', 3, ...
%         'NumSpaceTimeStreams', 3); 
%  
%     % Generate L-LTF and L-SIG field signals
%     txLLTF = wlanLLTF(cfgVHT); 
%     [txLSIG, LSIGBits] = wlanLSIG(cfgVHT);
% 
%     % Pass through a 3 x 2 quasi-static fading channel with AWGN 
%     H = 1/sqrt(2)*complex(randn(3, 2), randn(3, 2));
%     rxLLTF = awgn(txLLTF * H, 10);
%     rxLSIG = awgn(txLSIG * H, 10);
% 
%     % Perform channel estimation based on L-LTF
%     demodLLTF = wlanLLTFDemodulate(rxLLTF, chanBW, 1);
%     chanEst = wlanLLTFChannelEstimate(demodLLTF, chanBW);
% 
%     % Recover information bits in L-SIG
%     recLSIGBits = wlanLSIGRecover(rxLSIG, chanEst, 0.1, chanBW);
%
%     % Compare against original information bits
%     disp(isequal(LSIGBits, recLSIGBits));
% 
%   See also wlanRecoveryConfig, wlanLSIG, wlanLLTF, wlanLLTFDemodulate,
%   wlanLLTFChannelEstimate.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

narginchk(4, 5);
nargoutchk(0, 4);

% Calculate CPE if requested
if nargout>3
    calculateCPE = true;
else
    calculateCPE = false;
end

% Validate channel bandwidth input
coder.internal.errorIf(~ischar(chanBW) || ...
    ~any(strcmp(chanBW, {'CBW5','CBW10','CBW20','CBW40','CBW80','CBW160'})), ...
    'wlan:wlanLSIGRecover:InvalidChBandwidth');

% Get OFDM configuration
[cfgOFDM,dataInd,pilotInd] = wlanGetOFDMConfig(chanBW, 'Long', 'Legacy');
FFTLen = cfgOFDM.FFTLength;

% Validate L-SIG field signal input
validateattributes(rxLSIG, {'double'}, {'2d','finite','nrows',FFTLen*5/4}, ...
    'rxLSIG', 'L-SIG field signal'); 
numRx = size(rxLSIG, 2);

% Validate channel estimates
numSD = 48;
validateattributes(chanEst, {'double'}, {'3d','finite'}, ...
    'chanEst', 'channel estimates'); 

% Cross validate inputs
numST = numel([dataInd; pilotInd]); % Total number of subcarriers
coder.internal.errorIf(size(chanEst, 1) ~= numST, ...
    'wlan:wlanLSIGRecover:InvalidChanEst1D', numST);
coder.internal.errorIf(size(chanEst, 2) ~= 1, ...
    'wlan:wlanLSIGRecover:InvalidChanEst2D');
coder.internal.errorIf(size(chanEst, 3) ~= numRx, ...
    'wlan:wlanLSIGRecover:InvalidChanEst3D');

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

% ofdmOutData is [48*num20, 1, numRx]
[ofdmOutData, ofdmOutPilots] = wlanOFDMDemodulate(rxLSIG, cfgOFDM, symOffset); 

% Pilot phase tracking
if calculateCPE==true || strcmp(pilotPhaseTracking, 'PreEQ')
    % Get reference pilots, from IEEE Std 802.11-2012, Eqn 20-14
    z = 0; % No offset as first symbol with pilots
    refPilots = nonHTPilots(1, z, chanBW);

    % Estimate CPE and phase correct symbols
    cpe = commonPhaseErrorEstimate(ofdmOutPilots, chanEstPilots, ...
        refPilots);
    if strcmp(pilotPhaseTracking, 'PreEQ')
        ofdmOutData = commonPhaseErrorCorrect(ofdmOutData, cpe);
    end
    if calculateCPE==true
        varargout{1} = cpe.'; % Permute to Nsym-by-1
    end
end

% Merge num20 channel estimates and demodulated symbols together for the
% repeated subcarriers for data carrying subcarriers
ofdmDataOutOne20MHz = reshape(permute(reshape(ofdmOutData, numSD, [], numRx), ...
    [1 3 2]), numSD, 1, []); % [48, 1, num20*numRx]
chanEstDataOne20MHz = reshape(permute(reshape(chanEstData, numSD, [], numRx), ...
    [1 3 2]), numSD, 1, []); % [48, 1, num20*numRx]
% Perform equalization
[eqDataSym, csiData] = wlanEqualize(ofdmDataOutOne20MHz, ...
    chanEstDataOne20MHz, eqMethod, noiseVarEst);

% Constellation demapping
demodOut = wlanConstellationDemodulate(eqDataSym(:), 1, noiseVarEst); 
% Need the (:) above for codegen
demodOut = demodOut .* csiData;

% Deinterleaving
deintlvrOut = wlanBCCDeinterleave(demodOut, 'NON_HT', 48, 1);

% BCC decoding
bits = wlanBCCDecode(deintlvrOut, '1/2', 24);

% Parity check & rate check (the 4th bit must be 1)
failCheck = (mod(sum(bits(1:17)), 2) ~= bits(18)) || (bits(4) ~= 1);

end

% [EOF]