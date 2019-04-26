function [bits, failCRC, eqDataSym, varargout] = wlanHTSIGRecover(...
				rxHTSIG, chanEst, noiseVarEst, chanBW, varargin)
%WLANHTSIGRECOVER Recover information bits in HT-SIG field
% 
%   [BITS, FAILCRC] = wlanHTSIGRecover(RXHTSIG, CHANEST, NOISEVAREST,
%   CHANBW) recovers the information bits in the HT-SIG field and performs
%   CRC check.
% 
%   BITS is an int8 column vector of length 48 containing the recovered
%   information bits.
%
%   FAILCRC is true if BITS fails the CRC check. It is a logical scalar.
% 
%   RXHTSIG is the received time-domain HT-SIG field signal. It is a
%   Ns-by-Nr matrix of real or complex values, where Ns represents the
%   number of time-domain samples in the HT-SIG field and Nr represents the
%   number of receive antennas. Ns can be greater than the HT-SIG field
%   length; in this case additional samples at the end of RXHTSIG are not
%   used.
% 
%   CHANEST is the estimated channel at data and pilot subcarriers based on
%   the L-LTF. It is a real or complex array of size Nst-by-1-by-Nr, where
%   Nst represents the total number of occupied subcarriers.  The singleton
%   dimension corresponds to the single transmitted stream in the L-LTF
%   which includes the combined cyclic shifts if multiple transmit antennas
%   are used.
%   
%   NOISEVAREST is the noise variance estimate. It is a real, nonnegative
%   scalar.
%
%   CHANBW is a string describing the channel bandwidth and must be 'CBW20'
%   or 'CBW40'.
% 
%   [BITS, FAILCRC] = wlanHTSIGRecover(..., CFGREC) allows different
%   algorithm options for information bit recovery via the input CFGREC,
%   which is a <a href="matlab:help('wlanRecoveryConfig')">wlanRecoveryConfig</a> configuration object. When the CFGREC
%   input is not specified, the default property values of the
%   <a href="matlab:help('wlanRecoveryConfig')">wlanRecoveryConfig</a> object are adopted in the recovery.
%
%   [..., EQDATASYM, CPE] = wlanHTSIGRecover(...) also returns the
%   equalized subcarriers and common phase error.
%
%   EQDATASYM is a 48-by-2 complex matrix containing the equalized symbols
%   at data carrying subcarriers. There are 48 data carrying subcarriers in
%   each of the 2 OFDM symbols which constitute the HT-SIG field.
%
%   CPE is a column vector of length 2 containing the common phase error
%   between each of the 2 received and expected OFDM symbols.
%
%   Example: 
%   %  Recover the information bits in HT-SIG field via channel estimation 
%   %  on L-LTF over a 2 x 4 quasi-static fading channel
%
%     % Configure a HT configuration object 
%     chanBW = 'CBW40';
%     cfgHT = wlanHTConfig('ChannelBandwidth', 'CBW40', ...
%         'NumTransmitAntennas', 2, 'NumSpaceTimeStreams', 2); 
%  
%     % Generate L-LTF and HT-SIG field signals
%     txLLTF  = wlanLLTF(cfgHT); 
%     txHTSIG = wlanHTSIG(cfgHT);
% 
%     % Pass through a 2 x 4 quasi-static fading channel with AWGN 
%     H = 1/sqrt(2)*complex(randn(2, 4), randn(2, 4));
%     rxLLTF  = awgn(txLLTF  *H, 10);
%     rxHTSIG = awgn(txHTSIG *H, 10);
% 
%     % Perform channel estimation based on L-LTF
%     demodLLTF = wlanLLTFDemodulate(rxLLTF, chanBW, 1);
%     chanEst = wlanLLTFChannelEstimate(demodLLTF, chanBW);
% 
%     % Recover information bits in HT-SIG
%     [recHTSIGBits, failCRC] = wlanHTSIGRecover(rxHTSIG, chanEst, ...
%            0.1, chanBW);
%
%     % Check CRC validation
%     failCRC
% 
%   See also wlanRecoveryConfig, wlanHTSIG, wlanLLTF, wlanLLTFDemodulate,
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
    ~any(strcmp(chanBW, {'CBW20','CBW40'})), ...
    'wlan:wlanHTSIGRecover:InvalidChanBW');

% Get OFDM configuration
[cfgOFDM,dataInd,pilotInd] = wlanGetOFDMConfig(chanBW, 'Long', 'Legacy');

% Validate HT-SIG field signal input
validateattributes(rxHTSIG, {'double'}, ...
    {'2d','finite','nrows',cfgOFDM.FFTLength*5/2}, ...
    'rxHTSIG', 'HT-SIG field signal'); 
numRx = size(rxHTSIG, 2);

% Validate channel estimates
numSD = 48;
validateattributes(chanEst, {'double'}, {'3d','finite'}, ...
    'chanEst', 'channel estimates'); 

% Cross validate inputs
numST = numel([dataInd; pilotInd]); % Total number of occupied subcarriers
coder.internal.errorIf(size(chanEst, 1) ~= numST, ...
    'wlan:wlanHTSIGRecover:InvalidChanEst1D', numST);
coder.internal.errorIf(size(chanEst, 2) ~= 1, ...
    'wlan:wlanHTSIGRecover:InvalidChanEst2D');
coder.internal.errorIf(size(chanEst, 3) ~= numRx, ...
    'wlan:wlanHTSIGRecover:InvalidChanEst3D');

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
[ofdmOutData, ofdmOutPilots] = wlanOFDMDemodulate(rxHTSIG, cfgOFDM, symOffset);

% Pilot phase tracking
if calculateCPE==true || strcmp(pilotPhaseTracking, 'PreEQ')
    % Get reference pilots, from IEEE Std 802.11-2012, Eqn 20-17
    z = 1; % Offset by 1 to account for L-SIG pilot symbol
    refPilots = nonHTPilots(2, z, chanBW);
     
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
% repeated subcarriers for data subcarriers
ofdmOutDataOne20MHz = reshape(permute(reshape(ofdmOutData, numSD, [], 2, numRx), ...
    [1 3 4 2]), numSD, 2, []); % [numSD, 2, num20*numRx]
chanEstDataOne20MHz = reshape(permute(reshape(chanEstData, numSD, [], numRx), ...
    [1 3 2]), numSD, 1, []);   % [numSD, 1, num20*numRx]
% Perform equalization
[eqDataSym, csiData] = wlanEqualize(ofdmOutDataOne20MHz, ...
    chanEstDataOne20MHz, eqMethod, noiseVarEst);

% Constellation demapping
demodOut = wlanConstellationDemodulate(eqDataSym, 1, noiseVarEst, pi/2);
demodOut = demodOut .* repmat(csiData, 1, 2);

% Deinterleaving
deintlvrOut = zeros(size(demodOut));
deintlvrOut(:, 1) = wlanBCCDeinterleave(demodOut(:, 1), 'NON_HT', 48, 1);
deintlvrOut(:, 2) = wlanBCCDeinterleave(demodOut(:, 2), 'NON_HT', 48, 1);

% BCC decoding 
bits = wlanBCCDecode(deintlvrOut(:), '1/2'); 

% CRC detection
[~, failCRC] = wlanCRCDetect(bits(1:42));

end

% [EOF]