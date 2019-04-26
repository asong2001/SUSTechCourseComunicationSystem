function [bits, eqDataSym, varargout] = wlanNonHTDataRecover( ...
    rxNonHTData, chanEst, noiseVarEst, cfgNonHT, varargin)
%wlanNonHTDataRecover Recover information bits from Non-HT Data field signal
%
%   BITS = wlanNonHTDataRecover(RXNONHTDATA, CHANEST, NOISEVAREST,
%   CFGNONHT) recovers the information bits in the Non-HT Data field for a
%   Non-HT OFDM format transmission.
%
%   BITS is an int8 column vector of length 8*CFGNONHT.PSDULength
%   containing the recovered information bits.
%
%   RXNONHTDATA is the received time-domain NonHT-Data field signal. It is
%   a Ns-by-Nr matrix of real or complex values, where Ns represents the
%   number of time-domain samples in the NonHT-Data field and Nr represents
%   the number of receive antennas. Ns can be greater than the NonHT-Data
%   field length; in this case redundant samples at the end of RXNONHTDATA
%   are not used.
%
%   CHANEST is the estimated channel at data and pilot subcarriers based on
%   the L-LTF. It is a real or complex array of size Nst-by-1-by-Nr, where
%   Nst represents the total number of occupied subcarriers. The singleton
%   dimension corresponds to the single transmitted stream in the L-LTF
%   which includes the combined cyclic shifts if multiple transmit antennas
%   are used.
%
%   NOISEVAREST is the noise variance estimate. It is a real, nonnegative
%   scalar.
%
%   CFGNONHT is the format configuration object of type <a href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a> 
%   that specifies the Non-HT format parameters. Only OFDM modulation is
%   supported.
% 
%   BITS = wlanNonHTDataRecover(..., CFGREC) allows different algorithm
%   options for data recovery via the input CFGREC, which is a
%   <a href="matlab:help('wlanRecoveryConfig')">wlanRecoveryConfig</a> configuration object. When the CFGREC input is not 
%   specified, the default property values of the  <a href="matlab:help('wlanRecoveryConfig')">wlanRecoveryConfig</a> object
%   are adopted in the recovery.
%
%   [..., EQDATASYM, CPE] = wlanNonHTDataRecover(...) also returns the
%   equalized subcarriers and common phase error.
%
%   EQDATASYM is a complex 48-by-Nsym matrix containing the equalized
%   symbols at data carrying subcarriers. There are 48 data carrying
%   subcarriers in the Non-HT Data field. Nsym represents the number of
%   OFDM symbols in the NonHT-Data field.
%
%   CPE is a column vector of length Nsym containing the common phase error
%   between each received and expected OFDM symbol.
%   
%   Example: 
%   %  Recover a NONHT-Data field signal through a SISO AWGN channel
%   %  using ZF equalization.
%  
%     cfgNonHT = wlanNonHTConfig('PSDULength', 1024);  % NonHT OFDM 
%     txBits = randi([0 1], 8*cfgNonHT.PSDULength, 1); % PSDU bits
%     tx = wlanNonHTData(txBits, cfgNonHT);       % NonHT-Data field signal
% 
%     % Add AWGN, with noise variance of 1
%     rx = awgn(tx, 1, 1);
%
%     % Configure recovery object
%     cfgRec = wlanRecoveryConfig('EqualizationMethod', 'ZF'); 
%     % Recover PSDU bits 
%     rxBits = wlanNonHTDataRecover(rx, ones(52,1), 1, cfgNonHT, cfgRec);
%   
%     [numerr, ber] = biterr(rxBits, txBits); % Compare bits
%     disp(ber)
%     
%   See also wlanNonHTConfig, wlanRecoveryConfig, wlanLLTFChannelEstimate, 
%       wlanNonHTData. 

%   Copyright 2015 The MathWorks, Inc.

%#codegen

narginchk(4,5);
nargoutchk(0,3);

% Calculate CPE if requested
if nargout>2
    calculateCPE = true;
else
    calculateCPE = false;
end

% NonHT configuration input self-validation
validateattributes(cfgNonHT, {'wlanNonHTConfig'}, ...
    {'scalar'}, mfilename, 'format configuration object');
% Only applicable for OFDM and DUP-OFDM modulations
coder.internal.errorIf( ~strcmp(cfgNonHT.Modulation, 'OFDM'), ...
                        'wlan:wlanNonHTDataRecover:InvalidModulation');
s = validateConfig(cfgNonHT);

% Validate rxNonHTData
validateattributes(rxNonHTData, {'double'}, {'2d','finite'}, ...
    'rxHTData', 'Non-HT OFDM Data field signal'); 
% Validate chanEst
validateattributes(chanEst, {'double'}, {'3d','finite'}, ...
    'chanEst', 'channel estimates'); 
% Validate noiseVarEst
validateattributes(noiseVarEst, {'double'}, ...
    {'real','scalar','nonnegative','finite'}, ...
    'noiseVarEst', 'noise variance estimate'); 

% Optional recovery configuration input validation
if nargin == 5
    validateattributes(varargin{1}, {'wlanRecoveryConfig'}, ...
    {'scalar'}, mfilename, 'recovery configuration object');

    symOffset = varargin{1}.OFDMSymbolOffset;
    eqMethod  = varargin{1}.EqualizationMethod; 
    pilotPhaseTracking = varargin{1}.PilotPhaseTracking;
else % use defaults
    symOffset = 0.75;
    eqMethod  = 'MMSE'; 
    pilotPhaseTracking = 'PreEQ';
end

numRx = size(rxNonHTData, 2);

mcsTable = getRateTable(cfgNonHT);
pNcbpssi = mcsTable.NCBPS;

numOFDMSym = s.NumDataSymbols;

% Get OFDM configuration
[cfgOFDM,dataInd,pilotInd] = wlanGetOFDMConfig(cfgNonHT.ChannelBandwidth, ...
                                               'Long', 'Legacy');

% Cross validate inputs
numST = numel([dataInd; pilotInd]); % Total number of occupied subcarriers
coder.internal.errorIf(size(chanEst, 1) ~= numST, ...
    'wlan:wlanNonHTDataRecover:InvalidNHTChanEst1D', numST);
coder.internal.errorIf(size(chanEst, 2) ~= 1, ...
    'wlan:wlanNonHTDataRecover:InvalidNHTChanEst2D');
coder.internal.errorIf(size(chanEst, 3) ~= numRx, ...
    'wlan:wlanNonHTDataRecover:InvalidNHTChanEst3D');

% Extract data and pilot subcarriers from channel estimate
chanEstData = chanEst(dataInd,:,:);
chanEstPilots = chanEst(pilotInd,:,:);

% Cross-validation between inputs
minInputLen = numOFDMSym*(cfgOFDM.FFTLength+cfgOFDM.CyclicPrefixLength);
coder.internal.errorIf(size(rxNonHTData, 1) < minInputLen, ...
    'wlan:wlanNonHTDataRecover:ShortNHTDataInput', minInputLen);

% Processing 
% OFDM Demodulation 
[ofdmDemodData,ofdmDemodPilots] = wlanOFDMDemodulate( ...
    rxNonHTData(1:minInputLen, :), cfgOFDM, symOffset);

% Pilot phase tracking
if calculateCPE==true || strcmp(pilotPhaseTracking, 'PreEQ')
    % Get reference pilots, from IEEE Std 802.11-2012, Eqn 18-22
    z = 1; % Offset by 1 to account for L-SIG pilot symbol
    refPilots = nonHTPilots(numOFDMSym, z);
    
    % Estimate CPE and phase correct symbols
    cpe = commonPhaseErrorEstimate(ofdmDemodPilots, chanEstPilots, ...
        refPilots);
    if strcmp(pilotPhaseTracking, 'PreEQ')
        ofdmDemodData = commonPhaseErrorCorrect(ofdmDemodData, cpe);
    end
    if calculateCPE==true
        varargout{1} = cpe.'; % Permute to Nsym-by-1
    end
end

% Equalization
[eqDataSym, csiData] = wlanEqualize(ofdmDemodData, chanEstData, ...
                                    eqMethod, noiseVarEst);

% LLR Demodulation
qamDemodOut = wlanConstellationDemodulate(eqDataSym, mcsTable.NBPSCS, ...
                                          noiseVarEst);
% Apply bit-wise CSI
qamDemodOut = bsxfun(@times, ...
              reshape(qamDemodOut, mcsTable.NBPSCS, [], numOFDMSym), ...
              reshape(csiData, 1, [])); % [Nbpscs Nsd Nsym]

% Deinterleave per OFDM symbol
deintlvrOut = zeros(pNcbpssi*numOFDMSym, 1);
for symIdx = 1:numOFDMSym 
    deintlvrOut((symIdx-1)*pNcbpssi+(1:pNcbpssi), 1) = ...
        wlanBCCDeinterleave(reshape(qamDemodOut(:,:,symIdx), [], 1), ...
                    'NON_HT', pNcbpssi, mcsTable.NBPSCS);        
end

% Channel decoding
nhtDecBits = wlanBCCDecode(deintlvrOut, mcsTable.Rate);

% Combine bits and Descramble with derived initial State
iniY = nhtDecBits(1:7,1); % Take the first 7 bits for scrambler init
scramInit = [iniY(3)+iniY(7) iniY(2)+iniY(6) iniY(1)+iniY(5) ...
    iniY(4)+iniY(3)+iniY(7) iniY(3)+iniY(2)+iniY(6) ...
    iniY(2)+iniY(1)+iniY(5) iniY(1)+iniY(4)+iniY(3)+iniY(7)];
scramInit = mod(scramInit, 2);

% Remove the padded and tail bits to output only the data bits
descramDataOut = wlanScramble(nhtDecBits(1:(16+8*cfgNonHT.PSDULength),:),...
                              scramInit);

%   Remove the 16 service bits
bits = descramDataOut(17:end);   

end