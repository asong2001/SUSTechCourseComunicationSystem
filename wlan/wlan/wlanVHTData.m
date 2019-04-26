function y = wlanVHTData(PSDU,cfgVHT,varargin)
%wlanVHTData VHT Data field processing of the PSDU input
% 
%   Y = wlanVHTData(PSDU,CFGVHT) generates the VHT format Data field
%   time-domain waveform for the input PHY Service Data Unit (PSDU).
%
%   Y is the time-domain VHT Data field signal. It is a complex matrix of
%   size Ns-by-Nt, where Ns represents the number of time-domain samples
%   and Nt represents the number of transmit antennas.
%
%   PSDU is the PHY service data unit input to the PHY. For single-user,
%   PSDU can be a double or int8 typed binary column vector of length
%   CFGVHT.PSDULength*8. Alternatively, PSDU can be a row cell array with
%   length equal to number of users. The ith element in the cell array must
%   be a double or int8 typed binary column vector of length
%   CFGVHT.PSDULength(i)*8.
%
%   CFGVHT is the format configuration object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a> which
%   specifies the parameters for the VHT format.
%
%   Y = wlanVHTData(...,SCRAMINIT) optionally allows specification of the
%   scrambler initialization. When not specified, it defaults to a value of
%   93. When specified, it can be a double or int8-typed integer scalar or
%   1-by-Nu row vector between 1 and 127, inclusive, where Nu represents
%   the number of users. Alternatively, SCRAMINIT can be a double or
%   int8-typed binary 7-by-1 column vector or 7-by-Nu matrix, without any
%   all-zero column. If it is a scalar or column vector, it applies to all
%   users. Otherwise, each user can have its own scrambler initialization
%   as indicated by the corresponding column.
% 
%   Example: 
%   %  Generate signal for a 80MHz VHT Data field for single-user single
%   %  transmit antenna configuration.
% 
%      cfgVHT = wlanVHTConfig('ChannelBandwidth', 'CBW80', ...
%                            'NumTransmitAntennas', 1, ...
%                            'NumSpaceTimeStreams', 1, ...
%                            'MCS', 4);
%      inpPSDU = randi([0 1], cfgVHT.PSDULength*8, 1);    % PSDU in bits
%      y = wlanVHTData(inpPSDU,cfgVHT);
%
%   %  Generate signal for a 20MHz VHT Data field for two-user multiple
%   %  transmit antenna configuration.
% 
%      cfgVHT = wlanVHTConfig('ChannelBandwidth', 'CBW20', ...
%                            'NumUsers', 2, ...
%                            'GroupID', 2, ...
%                            'NumTransmitAntennas', 2, ...
%                            'NumSpaceTimeStreams', [1 1], ...
%                            'MCS', [4 8],...
%                            'APEPLength', [1024 2048]);
%      inp1PSDU = randi([0 1], cfgVHT.PSDULength(1)*8, 1); % User 1 PSDU 
%      inp2PSDU = randi([0 1], cfgVHT.PSDULength(2)*8, 1); % User 2 PSDU
%      inpPSDUs = {inp1PSDU, inp2PSDU};               % Concatenate PSDUs
%      y = wlanVHTData(inpPSDUs,cfgVHT);
%    
%   See also wlanVHTConfig, wlanWaveformGenerator, wlanVHTDataRecover.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

narginchk(2,3);
if nargin == 2
    scramInit = 93;     % As per IEEE Std 802.11-2012 Section L.1.5.2.
else
    scramInit = varargin{1};
end

% Validate format configuration input
validateattributes(cfgVHT, {'wlanVHTConfig'}, {'scalar'}, mfilename, ...
    'VHT format configuration object');
cfgInfo = validateConfig(cfgVHT, 'SMappingMCS');
numUsers = cfgVHT.NumUsers;

% Early return for NDP, if invoked by user.
if isscalar(cfgVHT.APEPLength) && (cfgVHT.APEPLength == 0)
    y = complex(zeros(0, cfgVHT.NumTransmitAntennas));
    return;
end

% Validate PSDU input
coder.internal.errorIf((numUsers > 1) && ~iscell(PSDU), ...
    'wlan:wlanVHTData:InvalidPSDUForMU');
if iscell(PSDU) % SU and MU
    validateattributes(PSDU, {'cell'}, {'row','numel',numUsers}, ...
        mfilename, 'PSDU input');
    
    for u = 1:numUsers
        validateattributes(PSDU{u}, {'double','int8'}, {'real', ...
            'integer','column','binary','numel',8*cfgVHT.PSDULength(u)}, ...
            mfilename, sprintf('PSDU input for user %d', u));
    end
    PSDUMU = PSDU;
else % SU
    validateattributes(PSDU, {'double','int8'}, {'real','integer', ...
        'column','binary','numel',8*cfgVHT.PSDULength(1)}, ...
        mfilename, 'PSDU input');
    PSDUMU = {PSDU};
end

% Validate scrambler init input
if isscalar(scramInit)      % [1 1]
    validateattributes(scramInit, {'double','int8'}, ...
        {'real','integer','nonempty','>=',1,'<=',127}, ...
        mfilename, 'Scrambler initialization');
    scramInitBits = repmat(de2bi(scramInit, 7, 'left-msb'), numUsers, 1);
elseif isrow(scramInit)     % [1 Nu]
    validateattributes(scramInit, {'double','int8'}, ...
        {'real','integer','nonempty','numel',numUsers,'>=',1,'<=',127}, ...
        mfilename, 'Scrambler initialization');
    scramInitBits = de2bi(scramInit, 7, 'left-msb');
elseif iscolumn(scramInit)  % [7, 1]
    validateattributes(scramInit, {'double','int8'}, ...
        {'real','integer','nonempty','numel',7,'binary'}, ...
        mfilename, 'Scrambler initialization');

    % Check for non-zero init
    coder.internal.errorIf(all(scramInit == 0), ...
        'wlan:wlanVHTData:InvalidScramInit');
    
    scramInitBits = repmat(scramInit', numUsers, 1);
else                        % [7 Nu]
    validateattributes(scramInit, {'double','int8'}, ...
        {'real','integer','nonempty','size',[7,numUsers],'binary'}, ...
        mfilename, 'Scrambler initialization');

    % Check for non-zero init
    coder.internal.errorIf(any(sum(scramInit) == 0), ...
        'wlan:wlanVHTData:InvalidScramInit');
    
    scramInitBits = scramInit';
end

% Set up implicit parameters
chanBW      = cfgVHT.ChannelBandwidth;
numSTSTotal = sum(cfgVHT.NumSpaceTimeStreams); 
numTx       = cfgVHT.NumTransmitAntennas;
numOFDMSym  = cfgInfo.NumDataSymbols; % [1 1]
numPadBits  = cfgInfo.NumPadBits;     % [1 Nu]
mcsTable    = getRateTable(cfgVHT);
numSD       = mcsTable.NSD(1);        % [1 1]
numSS       = mcsTable.Nss;           % [1 Nu]
vecAPEPLen  = repmat(cfgVHT.APEPLength, 1, ...
                     numUsers/length(cfgVHT.APEPLength)); % [1 Nu]
vecMCS      = repmat(cfgVHT.MCS, 1, numUsers/length(cfgVHT.MCS)); % [1 Nu]

% Get data subcarriers for all users
data = complex(zeros(numSD, numOFDMSym, sum(numSS)));
for u = 1:numUsers
    data(:, :, sum(numSS(1:u-1))+(1:numSS(u))) = ...        
        getDataSubcarrierPerUser(chanBW, vecAPEPLen(u), vecMCS(u), ...
        numOFDMSym, numPadBits(u), mcsTable, u, PSDUMU{u}, ...
        scramInitBits(u,:));
end

% STBC encoding, Section 22.3.10.9.4, IEEE Std 802.11ac-2013
if (numUsers == 1) && cfgVHT.STBC 
    stbcData = wlanSTBCEncode(data, numSTSTotal);
else
    stbcData = data;
end
    
% Get OFDM parameters
cfgOFDM = wlanGetOFDMConfig(chanBW, cfgVHT.GuardInterval, 'VHT', numSTSTotal);
FFTLen = cfgOFDM.FFTLength;
CPLen  = cfgOFDM.CyclicPrefixLength;

% Generate pilots for VHT, IEEE Std 802.11ac-2013, Eqn 22-95
z = 4; % Offset by 4 to allow for L-SIG, VHT-SIG-A, VHT-SIG-B pilot symbols
pilots = vhtPilots(numOFDMSym, z, chanBW, numSTSTotal);
% Permute to Nsp-by-Nsts-by-Nsym for efficient accessing
pilots = permute(pilots, [1 3 2]);

wout = complex(zeros((FFTLen + CPLen)*numOFDMSym, numTx));
packedData = complex(zeros(FFTLen, numSTSTotal));
csh = getCyclicShiftVal('VHT', numSTSTotal, FFTLen/64*20);

for i = 1:numOFDMSym
    % Data packing with pilot insertion
    packedData(cfgOFDM.DataIndices,  :) = stbcData(:,i,:);    
    packedData(cfgOFDM.PilotIndices, :) = pilots(:,:,i);
    
    % Tone rotation
    rotatedData = bsxfun(@times, packedData, cfgOFDM.CarrierRotations);
    
    % Cyclic shift
    dataCycShift = wlanCyclicShift(rotatedData, csh, FFTLen, 'Tx');
    
    % Spatial mapping
    dataSpatialMapped = wlanSpatialMapping(dataCycShift, ...
        cfgVHT.SpatialMapping, numTx, cfgVHT.SpatialMappingMatrix);
    
    % OFDM modulate - CPLen
    wout((1:(FFTLen+CPLen)).' + (i-1)*(FFTLen+CPLen), :) = ...
        wlanOFDMModulate(reshape(dataSpatialMapped, FFTLen, 1, numTx), ...
        CPLen);
end

% Scale and output
y = wout * cfgOFDM.NormalizationFactor;

end


function dataPerUser = getDataSubcarrierPerUser(chanBW, APEPLength, MCS, ...
    numOFDMSym, numPadBits, mcsTable, userIdx, PSDU, scramInitBits)
% Only for Data packets, not for Null-Data-Packets

chanBWMHz = str2double(chanBW(4:end));
numES     = mcsTable.NES(userIdx);
numSS     = mcsTable.Nss(userIdx);
rate      = mcsTable.Rate(userIdx);
numBPSCS  = mcsTable.NBPSCS(userIdx);
numCBPS   = mcsTable.NCBPS(userIdx);
numSeg    = 1 + strcmp(chanBW, 'CBW160');
numCBPSSI = numCBPS/numSS/numSeg;

% Determine VHT-SIG-B bits per user, 
%   Section 22.3.8.3.6, IEEE Std 802.11ac-2013
APEPLenOver4 = ceil(APEPLength/4);
if length(mcsTable.Nss) == 1 % SU PPDU allocation
    switch chanBW
      case 'CBW20' % 20
        bitsPerUser = [de2bi(APEPLenOver4, 17) ones(1, 3, 'int8')].';
      case 'CBW40' % 21
        bitsPerUser = [de2bi(APEPLenOver4, 19) ones(1, 2, 'int8')].';
      otherwise    % 23 for {'CBW80', 'CBW80+80', 'CBW160'}
        bitsPerUser = [de2bi(APEPLenOver4, 21) ones(1, 2, 'int8')].';
    end
else % MU PPDU allocation
    switch chanBW
      case 'CBW20' % 20
        bitsPerUser = int8([de2bi(APEPLenOver4, 16), de2bi(MCS, 4)].');
      case 'CBW40' % 21
        bitsPerUser = int8([de2bi(APEPLenOver4, 17), de2bi(MCS, 4)].');
      otherwise    % 23 for {'CBW80', 'CBW80+80', 'CBW160'} 
        bitsPerUser = int8([de2bi(APEPLenOver4, 19), de2bi(MCS, 4)].');
    end
end

% Assemble the service bits,
%   Section 22.3.10.2 and 22.3.10.3, IEEE Std 802.11ac-2013
%   SERVICE = [scrambler init = 0; Reserved = 0; CRC of SIGB bits];
serviceBits = [zeros(7,1); 0; wlanCRCGenerate(bitsPerUser)];

% Scramble padded data, Section 22.3.10.4, IEEE Std 802.11ac-2013
paddedData = [serviceBits; PSDU; zeros(numPadBits, 1)];
scrambData = wlanScramble(paddedData, scramInitBits);

% BCC Encoding
%   Reshape scrambled data as per IEEE Std 802.11ac-2013, Eq. 22-60
%   for multiple encoders
numTailBits = 6;
encodedStreams = reshape(scrambData, numES, []).';
encodedData = zeros(round((size(encodedStreams,1)+numTailBits)/rate),numES);
for i = 1:numES
    % BCC encoding, Section 22.3.10.5.3, IEEE Std 802.11ac-2013
    encodedData(:, i) = wlanBCCEncode([encodedStreams(:, i); ...
                                       zeros(numTailBits,1)], rate);
end

% Parse encoded data into streams
%   Section 22.3.10.6, IEEE Std 802.11ac-2013
streamParsedData = wlanStreamParser(encodedData, numSS, numBPSCS);

% Segment parsing for 160, 80+80 MHz
%   Section 22.3.10.7, IEEE Std 802.11ac-2013
parsedData = wlanSegmentParser(streamParsedData, chanBW, numBPSCS, ...
                               numCBPS, numES, 'Tx');      

% Interleaving, constellation mapping and segment deparsing
interleavedData = zeros(numCBPSSI, numSS, numSeg);
mappedData = complex(zeros(numCBPSSI/numBPSCS, numSS, numSeg));
%   [Nsd, Nsym, Nss]
dataPerUser = complex(zeros(numCBPSSI/numBPSCS*numSeg, numOFDMSym, numSS));
for i = 1:numOFDMSym    
    for segIdx = 1:numSeg % Frequency segmentation, if needed
        % Interleaving, Section 22.3.10.8, IEEE Std 802.11ac-2013
        interleavedData(:, :, segIdx) = wlanBCCInterleave( ...
            parsedData((i-1)*numCBPSSI + (1:numCBPSSI).', :, segIdx), ...
            'VHT', numCBPSSI, numBPSCS, chanBWMHz, numSS);
        
        % Constellation mapping, Section 22.3.10.9, IEEE Std 802.11ac-2013
        mappedData(:, :, segIdx) = wlanConstellationMapper( ...
            interleavedData(:, :, segIdx), numBPSCS);
    end
    
    % Frequency deparsing, Section 22.3.10.9.3, IEEE Std 802.11ac-2013
    dataPerUser(:,i,:) = reshape(wlanSegmentDeparser(mappedData, ...
                                 chanBW, 'Tx'), [], 1, numSS);
end

end

% [EOF]
