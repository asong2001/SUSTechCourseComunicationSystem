function y = wlanNonHTData(PSDU,cfgNonHT,varargin)
%WLANNONHTDATA Non-HT Data field processing of the PSDU
%
%   Y = wlanNonHTData(PSDU,CFGNONHT) generates the Non-HT format Data
%   field time-domain waveform for the input PLCP Service Data Unit (PSDU). 
%
%   Y is the time-domain Non-HT Data field signal. It is a complex matrix
%   of size Ns-by-Nt, where Ns represents the number of time-domain samples
%   and Nt represents the number of transmit antennas.
%
%   PSDU is the PHY service data unit input to the PHY. It is a double
%   or int8 typed column vector of length CFGNONHT.PSDULength*8, with each
%   element representing a bit.
%
%   CFGNONHT is the format configuration object of type <a href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a> which
%   specifies the parameters for the Non-HT format. Only OFDM modulation
%   type is supported.
%
%   Y = wlanNonHTData(...,SCRAMINIT) optionally allows specification of the
%   scrambler initialization, SCRAMINIT for the Data field. When not
%   specified, it defaults to a value of 93. When specified, it can be a
%   double or int8-typed positive scalar less than or equal to 127 or a
%   corresponding double or int8-typed binary vector of length 7.
%
%   Example:
%   %  Generate the signal for a 20 MHz Non-HT OFDM data field for 36 Mbps.
%
%     cfgNonHT = wlanNonHTConfig('MCS', 5);               % Configuration
%     inpPSDU = randi([0 1], cfgNonHT.PSDULength*8,1);    % PSDU in bits
%     y = wlanNonHTData(inpPSDU,cfgNonHT);
%   
%   See also wlanNonHTConfig, wlanLSIG, wlanNonHTDataRecover.

%   Copyright 2015 The MathWorks, Inc.

%#codegen
% For only NON-HT OFDM modulation.

narginchk(2,3);
if nargin==2
    scramInit = 93; % Default
else
    scramInit = varargin{1};
end

% Validate inputs
% Validate the format configuration object
validateattributes(cfgNonHT, {'wlanNonHTConfig'}, ...
    {'scalar'}, mfilename, 'format configuration object');
% Only applicable for OFDM and DUP-OFDM modulations
coder.internal.errorIf( ~strcmp(cfgNonHT.Modulation, 'OFDM'), ...
                        'wlan:wlanNonHTData:InvalidModulation');
s = validateConfig(cfgNonHT);
                    
if isscalar(scramInit)
    validateattributes(scramInit, {'double', 'int8'}, ...
    {'real','integer','scalar','>',0,'<=',127}, ...
    'wlanNonHTData:scramInit', 'Scrambler initialization');

    scramInitBits = de2bi(scramInit, 7, 'left-msb');
else
    validateattributes(scramInit, {'double', 'int8'}, ...
    {'real', 'integer', 'binary', 'vector', 'numel', 7}, ...
    'wlanNonHTData:scramInit', 'Scrambler initialization');

    % Check for non-zero, binary vector
    coder.internal.errorIf( all(scramInit==0), ...
        'wlan:wlanNonHTData:InvalidScramInit');

    scramInitBits = scramInit;
end

validateattributes(PSDU, {'double', 'int8'},...
    {'real', 'binary', 'size', [cfgNonHT.PSDULength*8 1]}, ...
    mfilename, 'PSDU input');

chanBW = cfgNonHT.ChannelBandwidth;

% Determine number of symbols and pad length
numSym = s.NumDataSymbols;
numPad = s.NumPadBits;
if (strcmp(chanBW, 'CBW10') || strcmp(chanBW, 'CBW5'))
    numTx = 1;  % override and set to 1 only, for 802.11j/p
else
    numTx  = cfgNonHT.NumTransmitAntennas;
end

mcsTable = getRateTable(cfgNonHT);
rate     = mcsTable.Rate;
numBPSCS = mcsTable.NBPSCS;
numCBPS  = mcsTable.NCBPS;
Ntail = 6;

cfgOFDM = wlanGetOFDMConfig(chanBW, 'Long', 'Legacy');
FFTLen = cfgOFDM.FFTLength;
CPLen  = cfgOFDM.CyclicPrefixLength;

%% Generate the data field

% SERVICE, Section 18.3.5.2, all zeros
serviceBits = zeros(16,1);

% Scramble padded data
%   [service; psdu; tail; pad] processing
paddedData = [serviceBits; PSDU; zeros(Ntail,1); zeros(numPad, 1)];
scrambData = wlanScramble(paddedData, scramInitBits);
% Zero-out the tail bits again for encoding
scrambData(16+length(PSDU) + (1:Ntail)) = zeros(Ntail,1);

% BCC Encoding
encodedData = wlanBCCEncode(scrambData, rate);

interleavedData = zeros(numCBPS, numSym);
for i = 1:numSym
    % Interleaving
    interleavedData(:,i) = wlanBCCInterleave(encodedData((i-1)*numCBPS + ...
                        (1:numCBPS).', :), 'NON_HT', numCBPS, numBPSCS);
end

% Constellation mapping
mappedData = wlanConstellationMapper(interleavedData, numBPSCS);

% Non-HT pilots, from IEEE Std 802.11-2012, Eqn 18-22
z = 1; % Offset by 1 to account for HT-SIG pilot symbol
pilotValues = nonHTPilots(numSym,z);
        
wout = complex(zeros((FFTLen + CPLen)*numSym, numTx));
packedData = complex(zeros(FFTLen, 1));
csh = getCyclicShiftVal('OFDM', numTx, 20);
for i = 1:numSym
    % Data packing with pilot insertion
    packedData(cfgOFDM.DataIndices, :) = mappedData(:,i);
        
    % Add pilots
    packedData(cfgOFDM.PilotIndices, :) = pilotValues(:,i);

    % Tone rotation and replicate over Tx
    pDataMat = repmat(packedData.*cfgOFDM.CarrierRotations, 1, numTx);
       
    % Cyclic shift applied per Tx
    dataCycShift =  wlanCyclicShift(pDataMat, csh, FFTLen, 'Tx');
        
    % OFDM modulate - pCPLen
    wout((1:(FFTLen+CPLen)).' + (i-1)*(FFTLen+CPLen), :) = ...
        wlanOFDMModulate(reshape(dataCycShift, FFTLen, 1, numTx), CPLen);
end

% Scale and output
y = wout * cfgOFDM.NormalizationFactor / sqrt(numTx);

end

% [EOF]
