function y = wlanHTData(PSDU,cfgHT,varargin)
%wlanHTData HT Data field processing of the PSDU input
%
%   Y = wlanHTData(PSDU, CFGHT) generates the HT-Mixed format Data field
%   time-domain waveform for the input Physical Layer Convergence Procedure
%   (PLCP) Service Data Unit (PSDU).
%
%   Y is the time-domain HT-Data field signal. It is a complex matrix of
%   size Ns-by-Nt, where Ns represents the number of time-domain samples
%   and Nt represents the number of transmit antennas.
%
%   PSDU is the PLCP service data unit input to the PHY. It is a double or
%   int8 typed column vector of length CFGHT.PSDULength*8, with each
%   element representing a bit.
%
%   CFGHT is the format configuration object of type <a href="matlab:help('wlanHTConfig')">wlanHTConfig</a> which
%   specifies the parameters for the HT-Mixed format.
%
%   Y = wlanHTData(..., SCRAMINIT) optionally allows specification of the
%   scrambler initialization for the Data field. When not specified, it
%   defaults to a value of 93. When specified, it can be a double or
%   int8-typed positive scalar less than or equal to 127 or a corresponding
%   double or int8-typed binary vector of length 7.
%
%   Example: Generate a signal for a MIMO 40MHz HT-Mixed data field.
%
%     cfgHT = wlanHTConfig('ChannelBandwidth', 'CBW40', ...
%                          'NumTransmitAntennas', 2, ...
%                          'NumSpaceTimeStreams', 2, ...
%                          'MCS', 12);
%     inpPSDU = randi([0 1], cfgHT.PSDULength*8, 1);    % PSDU in bits
%     y = wlanHTData(inpPSDU, cfgHT);
%   
%   See also wlanHTConfig, wlanWaveformGenerator, wlanHTDataRecover.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

narginchk(2,3);
if nargin==2
    scramInit = 93; % Default
else
    scramInit = varargin{1};
end

% Validate inputs
validateattributes(cfgHT, {'wlanHTConfig'}, {'scalar'}, mfilename, ...
                   'HT-Mixed format configuration object');
% Check dependent properties 
s = validateConfig(cfgHT);
               
if isscalar(scramInit)
    validateattributes(scramInit, {'double', 'int8'}, ...
    {'real', 'integer', 'scalar', '>', 0, '<=', 127}, ...
    'wlanHTData:scramInit', 'Scrambler initialization');

    scramInitBits = de2bi(scramInit, 7, 'left-msb');
else
    validateattributes(scramInit, {'double', 'int8'}, ...
    {'real', 'integer', 'binary', 'vector', 'numel', 7}, ...
    'wlanHTData:scramInit', 'Scrambler initialization');

    % Check for non-zero binary vector
    coder.internal.errorIf( all(scramInit==0), ...
        'wlan:wlanHTData:InvalidScramInit');

    scramInitBits = scramInit;
end

% Validate PSDU input
validateattributes(PSDU, {'double', 'int8'}, ...
    {'real', 'integer', 'binary', 'column', 'size', [cfgHT.PSDULength*8 1]}, ...
    mfilename, 'PSDU input');

numTx   = cfgHT.NumTransmitAntennas;
if cfgHT.PSDULength == 0
    y = complex(zeros(0, numTx));
    return;
end

% Get number of symbols and pad length
numSym = s.NumDataSymbols;
numPad = s.NumPadBits;

mcsTable = getRateTable(cfgHT);
numNes   = mcsTable.NES;
rate     = mcsTable.Rate;
numBPSCS = mcsTable.NBPSCS;
numCBPS  = mcsTable.NCBPS;
numSS    = mcsTable.Nss;
numSTS   = cfgHT.NumSpaceTimeStreams;
Ntail    = 6;
numCBPSSI = numCBPS/numSS;

% Get OFDM parameters
cfgOFDM = wlanGetOFDMConfig(cfgHT.ChannelBandwidth, ...
    cfgHT.GuardInterval, 'HT', numSTS);
FFTLen = cfgOFDM.FFTLength;
CPLen  = cfgOFDM.CyclicPrefixLength;
chanBWInMHz = FFTLen/64 * 20;

%%
% Generate the data field

% SERVICE bits, all zeros, IEEE Std 802.11-2012, Section 20.3.11.2
serviceBits = zeros(16,1);

% Scramble padded data
%   [service; psdu; tail; pad] processing
paddedData = [serviceBits; PSDU; zeros(Ntail*numNes,1); zeros(numPad, 1)];
scrambData = wlanScramble(paddedData, scramInitBits);
% Zero-out the tail bits again for encoding
scrambData(16+length(PSDU) + (1:Ntail*numNes)) = zeros(Ntail*numNes,1);

% BCC Encoding
%   Reshape scrambled data as per IEEE Std 802.11-2012, Eq 20-33, 
%   for multiple encoders
encodedStreams = reshape(scrambData, numNes, []).';
encodedData = zeros(round(size(encodedStreams,1)/rate), numNes, 'int8');
for i = 1:numNes
    encodedData(:, i) = wlanBCCEncode(encodedStreams(: , i), rate);
end

% Parse encoded data into streams
streamParsedData = wlanStreamParser(encodedData, numSS, numBPSCS);

wout = complex(zeros((FFTLen + CPLen)*numSym, numTx));
packedData = complex(zeros(FFTLen, numSTS));
mappedData = complex(zeros(numCBPSSI/numBPSCS, numSym, numSS));
csh = getCyclicShiftVal('VHT', numSTS, chanBWInMHz);

for i = 1:numSym
    % Interleaving
    interleavedData = wlanBCCInterleave( ...
        streamParsedData((i-1)*numCBPSSI + (1:numCBPSSI).', :), ...
        'HT', numCBPSSI, numBPSCS, chanBWInMHz, numSS);

    % Constellation mapping
    mappedData(:,i,:) = wlanConstellationMapper(interleavedData, numBPSCS);
end

if numSTS > numSS
    stbcData = wlanSTBCEncode(mappedData, numSTS);
else
    stbcData = mappedData;
end

% Generate pilots for HT, IEEE Std 802.11-2012, Eqn 22-58/59
z = 3; % offset by 3 to allow for L-SIG, HT-SIG pilot symbols, Eqn 20-58
pilots = htPilots(numSym,z,cfgHT.ChannelBandwidth,numSTS);
% Permute to Nsp-by-Nsts-by-Nsym for efficient accessing
pilots = permute(pilots,[1 3 2]);

for i = 1:numSym
    % Data packing with pilot insertion
    packedData(cfgOFDM.DataIndices,  :) = stbcData(:,i,:);
    packedData(cfgOFDM.PilotIndices, :) = pilots(:,:,i);

    % Tone rotation
    rotatedData = bsxfun(@times, packedData, cfgOFDM.CarrierRotations);
    
    % Cyclic shift applied per STS
    dataCycShift =  wlanCyclicShift(rotatedData, csh, FFTLen, 'Tx');
    
    % Spatial mapping
    dataSpMapped = wlanSpatialMapping(dataCycShift,... 
        cfgHT.SpatialMapping, numTx, cfgHT.SpatialMappingMatrix);
    
    % OFDM modulate
    wout((1:(FFTLen+CPLen)).' + (i-1)*(FFTLen+CPLen), :) = ...
        wlanOFDMModulate(reshape(dataSpMapped, FFTLen, 1, numTx), CPLen);
end

% Scale and output
y = wout * cfgOFDM.NormalizationFactor;

end

% [EOF]
