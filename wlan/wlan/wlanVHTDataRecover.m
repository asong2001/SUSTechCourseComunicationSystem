function [bits, CRCBits, eqDataSym, varargout] = wlanVHTDataRecover( ...
    rxVHTData, chanEst, noiseVarEst, cfgVHT, varargin)
%WLANVHTDATARECOVER Recover bits from VHT Data field signal
% 
%   [BITS, CRCBITS] = wlanVHTDataRecover(RXVHTDATA, CHANEST, NOISEVAREST,
%   CFGVHT) recovers the bits in the VHT-Data field for a VHT format
%   transmission.
%
%   BITS is an int8 column vector of length 8*CFGVHT.PSDULength containing
%   the recovered information bits.
%
%   CRCBITS is an int8 column vector of length 8 containing the VHT-Data
%   field checksum bits.
%
%   RXVHTDATA is the received time-domain VHT Data field signal. It is a
%   Ns-by-Nr matrix of real or complex values, where Ns represents the
%   number of time-domain samples in the VHT Data field and Nr represents
%   the number of receive antennas. Ns can be greater than the VHT Data
%   field length; in this case additional samples at the end of RXVHTDATA
%   are not used.
% 
%   CHANEST is the estimated channel at data and pilot subcarriers based on
%   the VHT-LTF. It is a real or complex array of size Nst-by-Nsts-by-Nr,
%   where Nst represents the total number of occupied subcarriers.
%
%   NOISEVAREST is the noise variance estimate. It is real, nonnegative
%   scalar.
%
%   CFGVHT is the format configuration object of type <a 
%   href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a>, which
%   specifies the parameters for the VHT format.
% 
%   [BITS,CRCBITS] = wlanVHTDataRecover(..., CFGREC) allows different
%   algorithm options for data recovery via the input CFGREC, which is a
%   <a href="matlab:help('wlanRecoveryConfig')">wlanRecoveryConfig</a> configuration object. When the CFGREC input is not
%   specified, the default property values of the <a href="matlab:help('wlanRecoveryConfig')">wlanRecoveryConfig</a> object
%   are adopted in the recovery.
%
%   [..., EQDATASYM, CPE] = wlanVHTDataRecover(...) also returns the
%   equalized subcarriers and common phase error.
%
%   EQDATASYM is a complex Nsd-by-Nsym-by-Nss array containing the
%   equalized symbols at data carrying subcarriers. Nsd represents the
%   number of data subcarriers, Nsym represents the number of OFDM symbols
%   in the VHT-Data field, and Nss represents the number of spatial
%   streams.
%
%   CPE is a column vector of length Nsym containing the common phase error
%   between each received and expected OFDM symbol.
%
%   Example: 
%   %  Recover bits in VHT Data field via channel estimation on VHT-LTF 
%   %  over a 2 x 2 quasi-static fading channel
%
%     % Configure a VHT configuration object 
%     chanBW = 'CBW160';
%     cfgVHT = wlanVHTConfig('ChannelBandwidth',    chanBW, ...
%         'NumTransmitAntennas', 2, 'NumSpaceTimeStreams', 2, ...
%         'APEPLength',          512); 
%  
%     % Generate VHT-LTF and VHT Data field signals
%     txDataBits = randi([0 1], 8*cfgVHT.PSDULength, 1);
%     txVHTLTF  = wlanVHTLTF(cfgVHT); 
%     txVHTData = wlanVHTData(txDataBits, cfgVHT);
% 
%     % Pass through a 2 x 2 quasi-static fading channel with AWGN 
%     H = 1/sqrt(2)*complex(randn(2, 2), randn(2, 2));
%     rxVHTLTF  = awgn(txVHTLTF  * H, 10);
%     rxVHTData = awgn(txVHTData * H, 10);
% 
%     % Perform channel estimation based on VHT-LTF
%     demodVHTLTF = wlanVHTLTFDemodulate(rxVHTLTF, cfgVHT, 1);
%     chanEst = wlanVHTLTFChannelEstimate(demodVHTLTF, cfgVHT);
% 
%     % Configure a recovery object using ZF equalization
%     cfgRec = wlanRecoveryConfig('EqualizationMethod', 'ZF'); 
% 
%     % Recover information bits in VHT Data
%     rxDataBits = wlanVHTDataRecover(rxVHTData, chanEst, 0.1, ...
%         cfgVHT, cfgRec);
%
%     % Compare against original information bits
%     disp(isequal(txDataBits, rxDataBits));
%
%   See also wlanVHTConfig, wlanRecoveryConfig, wlanVHTData, wlanVHTLTF,
%   wlanVHTLTFDemodulate, wlanVHTLTFChannelEstimate.

%   Copyright 2015 The MathWorks, Inc.

%#codegen
%#ok<*EMCA>

narginchk(4,5);
nargoutchk(0,4);

% Calculate CPE if requested
if nargout>3
    calculateCPE = true;
else
    calculateCPE = false;
end

% VHT configuration input self-validation
validateattributes(cfgVHT, {'wlanVHTConfig'}, {'scalar'}, ...
    mfilename, 'VHT format configuration object');
coder.internal.errorIf(cfgVHT.NumUsers ~= 1, ...
    'wlan:wlanVHTDataRecover:SingleUserOnly');
cfgInfo = validateConfig(cfgVHT, 'MCS');
mcsTable = getRateTable(cfgVHT);
chanBW = cfgVHT.ChannelBandwidth;

if cfgVHT.PSDULength(1) == 0
    bits     = zeros(0, 1, 'int8');
    CRCBits  = zeros(0, 1, 'int8');
    eqDataSym   = zeros(mcsTable.NSD, 0, mcsTable.Nss);
    if calculateCPE==true
        varargout{1} = []; % CPE
    end
    return;
end

% Optional recovery configuration input validation
if nargin == 5
    validateattributes(varargin{1}, {'wlanRecoveryConfig'}, {'scalar'}, ...
    mfilename, 'recovery configuration object');

    symOffset = varargin{1}.OFDMSymbolOffset;
    eqMethod  = varargin{1}.EqualizationMethod;
    pilotPhaseTracking  = varargin{1}.PilotPhaseTracking; 
else
    symOffset = 0.75;
    eqMethod  = 'MMSE';
    pilotPhaseTracking = 'PreEQ';
end

% Signal input self-validation
validateattributes(rxVHTData,   {'double'}, {'2d','finite'}, ...
    'rxVHTData',   'VHT-Data field signal'); 
validateattributes(chanEst, {'double'}, {'3d','finite'}, ...
    'chanEst', 'channel estimation'); 
validateattributes(noiseVarEst, {'double'}, ...
    {'real','scalar','nonnegative','finite'}, ...
    'noiseVarEst', 'noise variance estimation'); 

% Set up some implicit configuration parameters
numBPSCS   = mcsTable.NBPSCS(1);    % Number of coded bits per single carrier
numCBPS    = mcsTable.NCBPS(1);     % Number of coded bits per OFDM symbol
numES      = mcsTable.NES(1);       % Number of encoded streams
numSS      = mcsTable.Nss(1);       % Number of spatial streams
numSTS     = cfgVHT.NumSpaceTimeStreams;
numSeg     = strcmp(chanBW, 'CBW160') + 1;
% Number of coded bits per OFDM symbol, per spatial stream, per segment
numCBPSSI  = numCBPS/numSS/numSeg;  
numRx      = size(rxVHTData, 2);
numOFDMSym = cfgInfo.NumDataSymbols(1);

% Get OFDM configuration
[cfgOFDM, dataInd, pilotInd] = wlanGetOFDMConfig(chanBW, ...
    cfgVHT.GuardInterval, 'VHT', sum(numSTS));

% Cross-validation between inputs
numST = numel([dataInd; pilotInd]); % Total number of occupied subcarriers
coder.internal.errorIf(size(chanEst, 1) ~= numST, ...
    'wlan:wlanVHTDataRecover:InvalidChanEst1D', numST);
coder.internal.errorIf(size(chanEst, 2) ~= numSTS, ...
    'wlan:wlanVHTDataRecover:InvalidChanEst2D', numSTS);
coder.internal.errorIf(size(chanEst, 3) ~= numRx, ...
    'wlan:wlanVHTDataRecover:InvalidChanEst3D');

% Extract data and pilot subcarriers from channel estimate
chanEstData = chanEst(dataInd,:,:);
chanEstPilots = chanEst(pilotInd,:,:);

% Cross-validation between inputs
minInputLen = numOFDMSym*(cfgOFDM.FFTLength+cfgOFDM.CyclicPrefixLength);
coder.internal.errorIf(size(rxVHTData, 1) < minInputLen, ...
    'wlan:wlanVHTDataRecover:ShortDataInput', minInputLen);

% OFDM demodulation
[ofdmDemodData, ofdmDemodPilots] = wlanOFDMDemodulate( ...
    rxVHTData(1:minInputLen, :), cfgOFDM, symOffset);

% Pilot phase tracking
if calculateCPE==true || strcmp(pilotPhaseTracking, 'PreEQ')
    % Get reference pilots, from Eqn 22-95, IEEE Std 802.11ac-2013
    % Offset by 4 to allow for L-SIG, VHT-SIG-A, VHT-SIG-B pilot symbols
    z = 4; 
    refPilots = vhtPilots(numOFDMSym, z, chanBW, numSTS);
    
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
if cfgVHT.STBC
    [eqDataSym, dataCSI] = wlanSTBCCombine(ofdmDemodData, chanEstData, ...
        numSS, eqMethod, noiseVarEst);
else
    [eqDataSym, dataCSI] = wlanEqualize(ofdmDemodData, ...
        chanEstData, eqMethod, noiseVarEst);
end

% Segment deparsing
%   [Nsd/Nseg Nsym Nss Nseg]
parserOut    = wlanSegmentDeparser(eqDataSym, chanBW, 'Rx'); 
%   [Nsd/Nseg 1 Nss Nseg]
csiParserOut = wlanSegmentDeparser(reshape(dataCSI, [], 1, numSS), ...
                                   chanBW, 'Rx');

deintlvrOut = zeros(numCBPSSI*numOFDMSym, numSS, numSeg);
chBW = str2double(chanBW(4:end)); % For deinterleaving

% Constellation demapping and deinterleaving per segment
for segIdx = 1 : numSeg 
    qamDemodOut = wlanConstellationDemodulate(parserOut(:,:,:,segIdx), ...
                      numBPSCS, noiseVarEst); % [Nbpscs*Nsd/Nseg Nsym Nss]    
    qamDemodOut = bsxfun(@times, ...
        reshape(qamDemodOut, numBPSCS, [], numOFDMSym, numSS), ...
        reshape(csiParserOut(:,:,:,segIdx), 1, [], 1, numSS)); 
        % [Nbpscs Nsd/Nseg Nsym Nss]
    
    for symIdx = 1:numOFDMSym
        deintlvrOut((symIdx-1)*numCBPSSI+(1:numCBPSSI), :, segIdx) = ...
            wlanBCCDeinterleave( ...
            reshape(real(qamDemodOut(:,:,symIdx,:)), [], numSS), ...
            'VHT', numCBPSSI, numBPSCS, chBW, numSS);
    end
end

% Segment parsing
segDeparserOut = wlanSegmentParser(deintlvrOut, chanBW, ...
                                   numBPSCS, numCBPS, numES, 'Rx'); 

% Stream deparsing
streamDeparserOut = wlanStreamDeparser(segDeparserOut, numES, numBPSCS);

% Channel decoding
chanDecOut = zeros(round(size(streamDeparserOut, 1)*mcsTable.Rate(1)), ...
    numES, 'int8');
for nesIdx = 1 : numES 
    chanDecOut(:, nesIdx) = wlanBCCDecode(streamDeparserOut(:, nesIdx), ...
                                          mcsTable.Rate(1));
end

% Descrambling with derived initial state
numTailBits = 6; 
preDescramBits = reshape(chanDecOut(1:end-numTailBits,:)', [], 1);
scramInit = mod([preDescramBits(3:-1:1) + preDescramBits(7:-1:5); ...
                 preDescramBits(4:-1:2) + preDescramBits(3:-1:1) + ...
                 preDescramBits(7:-1:5); sum(preDescramBits([1 3 4 7]))]', 2);
descramBits = wlanScramble(preDescramBits(1:16+8*cfgVHT.PSDULength(1)), ...
    scramInit);

% Outputs
CRCBits = descramBits(9:16);
bits = descramBits(17:end);

end

% [EOF]