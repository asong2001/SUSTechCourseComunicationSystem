function txWaveform = wlanWaveformGenerator(dataBits,cfgFormat,varargin)
% wlanWaveformGenerator WLAN waveform generation
%   WAVEFORM = wlanWaveformGenerator(DATA,CFGFORMAT) generates a waveform
%   for a given format configuration and information bits. The generated
%   waveform contains a single packet with no idle time. For OFDM based
%   formats, the data scrambler initial states is 93 and the packet is
%   windowed for spectral controls with a windowing transition time of 1e-7
%   seconds.
%
%   WAVEFORM is a complex Ns-by-Nt matrix containing the generated
%   waveform, where Ns is the number of time domain samples, and Nt is the
%   number of transmit antennas.
%
%   DATA is the information bits including any MAC padding to be coded
%   across the number of packets to generate, i.e., representing multiple
%   concatenated PSDUs. It can be a double or int8 typed binary vector.
%   Alternatively, it can be a scalar cell array or a vector cell array
%   with length equal to number of users. Each element of the cell array
%   must be a double or int8 typed, binary vector. When DATA is a vector or
%   scalar cell array, it applies to all users. When DATA is a vector cell
%   array, each element applies to a single user. For each user, the bit
%   vector applied is looped if the number of bits required across all
%   packets of the generation exceeds the length of the vector provided.
%   This allows a short pattern to be entered, e.g. [1;0;0;1]. This pattern
%   will be repeated as the input to the PSDU coding across packets and
%   users. The number of data bits taken from a data stream for the ith
%   user when generating a packet is given by the ith element of the
%   CFGFORMAT.PSDULength property times eight.
%
%   CFGFORMAT is a format configuration object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a>, 
%   <a href="matlab:help('wlanHTConfig')">wlanHTConfig</a>, or <a href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a>. The format of the generated waveform 
%   is determined by the type of CFGFORMAT. The properties of CFGFORMAT are
%   used to parameterize the packets generated including the data rate and
%   PSDU length.
% 
%   WAVEFORM = wlanWaveformGenerator(DATA,CFGFORMAT,Name,Value) specifies
%   additional name-value pair arguments described below. When a name-value
%   pair is not specified, its default value is used.
%   
%   'NumPackets'               The number of packets to generate. It must
%                              be a positive integer. The default value is
%                              1.
%   
%   'IdleTime'                 The length in seconds of an idle period
%                              after each generated packet. It must be 0 or
%                              greater than or equal to 2e-6 seconds. The
%                              default value is 0 seconds.
% 
%   'ScramblerInitialization'  Scrambler initial state(s), applied for OFDM
%                              based formats. It must be a double or
%                              int8-typed scalar or matrix containing
%                              integer values between 1 and 127 inclusive.
%                              If a scalar is provided all packets are
%                              initialized with the same state for all
%                              users. Specifying a matrix allows a
%                              different initial state to be used per user
%                              and per packet. Each column specifies the
%                              initial states for a single user. If a
%                              single column is provided, the same initial
%                              states will be used for all users. Each row
%                              represents the initial state of each packet
%                              to generate. Internally the rows are looped
%                              if the number of packets to generate exceeds
%                              the number of rows of the matrix provided.
%                              The default value is 93, which is the
%                              example state given in IEEE Std 802.11-2012
%                              Section L.1.5.2.
% 
%   'WindowTransitionTime'     The windowing transition length in seconds,
%                              applied to OFDM based formats. It must be a
%                              nonnegative scalar and no greater than
%                              6.4e-6 seconds. Specifying it as 0 turns off
%                              windowing. The default value is 1e-7
%                              seconds.
%   Examples:
%
%   Example 1:
%       %  Generate a time domain signal txWaveform for an 802.11ac VHT
%       %  transmission with 10 packets and 20 microsecond idle period 
%       %  between packets.
%
%       numPkts = 10;                   % 10 packets in the waveform
%      
%       cfgVHT = wlanVHTConfig();       % Create format configuration
%       % Change properties from defaults
%       cfgVHT.NumTransmitAntennas = 2; % 2 transmit antennas
%       cfgVHT.NumSpaceTimeStreams = 2; % 2 spatial streams
%       cfgVHT.MCS = 1;                 % Modulation: QPSK Rate: 1/2 
%       cfgVHT.APEPLength = 1024;       % A-MPDU length in bytes
%
%       % Create bit vector containing concatenated PSDUs
%       numBits = cfgVHT.PSDULength*8*numPkts;
%       dataBits = randi([0 1],numBits,1);
%
%       txWaveform = wlanWaveformGenerator(dataBits, cfgVHT, ...
%           'NumPackets', numPkts, 'IdleTime', 20e-6, ...
%           'WindowTransitionTime', 1e-7);
%
%   Example 2:
%       %  Produce a waveform containing a single 802.11a packet without 
%       %  any windowing.
%          
%       cfgNonHT = wlanNonHTConfig(); % Create format configuration
%
%       psdu = randi([0 1], cfgNonHT.PSDULength*8, 1); % Create a PSDU
%
%       txWaveform = wlanWaveformGenerator(psdu, cfgNonHT, ...
%           'WindowTransitionTime', 0);  % Disable windowing
%
%   See also wlanVHTConfig, wlanHTConfig, wlanNonHTConfig.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

% Check number of input arguments
coder.internal.errorIf((nargin ~= 3) && mod(nargin, 2) == 1, ...
    'wlan:wlanWaveformGenerator:InvalidNumInputs');

% Validate the format configuration object is a valid type
validateattributes(cfgFormat, ...
    {'wlanVHTConfig','wlanHTConfig','wlanNonHTConfig'}, ...
    {'scalar'}, mfilename, 'format configuration object');
s = validateConfig(cfgFormat);
inDSSSMode = isa(cfgFormat,'wlanNonHTConfig') && ...
    strcmpi(cfgFormat.Modulation,'DSSS');

% Set maximum limits for windowing transition time based on bandwidth and
% format for the validity check performed later.
maxWinTransitTime = 6.4e-6;         % default
maxWinTransitTimeStr = '6.4e-06';   % string for codegen
if isa(cfgFormat, 'wlanNonHTConfig')
    if strcmpi(cfgFormat.Modulation, 'OFDM')
        switch cfgFormat.ChannelBandwidth
            case 'CBW5'
                maxWinTransitTime = 6.4e-6; % seconds
                maxWinTransitTimeStr = '6.4e-06'; 
            case 'CBW10'
                maxWinTransitTime = 3.2e-6; % seconds
                maxWinTransitTimeStr = '3.2e-06'; 
            otherwise  % 'CBW20'
                maxWinTransitTime = 1.6e-6; % seconds
                maxWinTransitTimeStr = '1.6e-06'; 
        end
    end
else % HT/VHT
    maxWinTransitTime = 1.6e-6; % seconds
    maxWinTransitTimeStr = '1.6e-06'; % default
end

if nargin == 3 % wlanGeneratorConfig
    cfgWaveGen = varargin{1};
    % Validate the waveform generator config object is the correct type
    validateattributes(cfgWaveGen, {'wlanGeneratorConfig'}, {'scalar'}, ...
       mfilename, 'waveform generator configuration object');
    numPackets     = cfgWaveGen.NumPackets;
    idleTime       = cfgWaveGen.IdleTime;
    scramblerInit  = cfgWaveGen.ScramblerInitialization;
    winTransitTime = cfgWaveGen.WindowTransitionTime;
    windowing      = (~inDSSSMode) && (cfgWaveGen.Windowing) && (winTransitTime > 0);
    coder.internal.errorIf(windowing && (winTransitTime > maxWinTransitTime), ...
                'wlan:wlanWaveformGenerator:ExceededMaxTransitionTime', ...
                maxWinTransitTimeStr);
else % P-V pairs
    % Define default P-V pair values
    numPackets     = 1;
    idleTime       = 0;
    scramblerInit  = 93;
    winTransitTime = 1e-7;
    
    % Validate each P-V pair
    for i = 3:2:nargin
        prop = varargin{i-2};
        val  = varargin{i-1};
        
        coder.internal.errorIf(~ischar(prop) || ...
            ~any(strcmp(prop, {'NumPackets', 'IdleTime', ...
            'ScramblerInitialization', 'WindowTransitionTime'})), ...
            'wlan:wlanWaveformGenerator:InvalidProperty');
        
        switch prop
          case 'NumPackets'
            validateattributes(val,{'numeric'},{'scalar','integer','>=',0}, ...
                mfilename,'''NumPackets'' value');
            numPackets = val;
          case 'IdleTime'
            validateattributes(val,{'numeric'}, ...
                {'scalar','real','>=',0},mfilename,'''IdleTime'' value');
            coder.internal.errorIf((val > 0) && (val < 2e-6), ...
                'wlan:wlanWaveformGenerator:InvalidIdleTimeValue');
            idleTime = val;
          case 'ScramblerInitialization'
            if ~inDSSSMode
                validateattributes(val,{'double','int8'}, ...
                    {'real','integer','2d','nonempty','>=',1,'<=',127}, ...
                    mfilename,'''ScramblerInitialization'' value');
            end
            scramblerInit = val;
          otherwise % 'WindowTransitionTime'
            if ~inDSSSMode
                validateattributes(val, {'numeric'}, ...
                    {'real','scalar','>=',0,'<=',maxWinTransitTime}, ...
                     mfilename,'''WindowTransitionTime'' value');
            end
            winTransitTime = val;
        end
    end
    windowing = (~inDSSSMode) && (winTransitTime > 0);
end

if isa(cfgFormat, 'wlanVHTConfig')
    numUsers = cfgFormat.NumUsers;
else
    numUsers = 1;
end

% Cross validation
coder.internal.errorIf( ...
    all(size(scramblerInit, 2) ~= [1 numUsers]), ...
    'wlan:wlanWaveformGenerator:ScramInitNotMatchNumUsers');

% Validate that data bits are present if PSDULength is nonzero
if iscell(dataBits) % SU and MU
    % Data must be a scalar cell or a vector cell of length Nu
    coder.internal.errorIf(~isvector(dataBits) || ...
        all(length(dataBits) ~= [1 numUsers]), ...
        'wlan:wlanWaveformGenerator:InvalidDataCell');
    
    for u = 1:length(dataBits)
        if ~isempty(dataBits{u}) && any(cfgFormat.PSDULength > 0) % Data packet
            validateattributes(dataBits{u}, {'double','int8'}, ...
                {'real','integer','vector','binary'}, ...
                mfilename, 'each element in cell data input');
        else
            % Empty data check if not NDP
            coder.internal.errorIf( ...
                any(cfgFormat.PSDULength > 0) && isempty(dataBits{u}), ...
                'wlan:wlanWaveformGenerator:NoData');
        end
    end
    if (numUsers > 1) && isscalar(dataBits) 
        % Columnize and expand to a [1 Nu] cell
        dataCell = repmat({int8(dataBits{1}(:))}, 1, numUsers);
    else % Columnize each element
        dataCell = repmat({int8(1)}, 1, numUsers); 
        for u = 1:numUsers                
            dataCell{u} = int8(dataBits{u}(:));
        end
    end
else % SU and MU: Data must be a vector
    if ~isempty(dataBits) && any(cfgFormat.PSDULength > 0) % Data packet
        validateattributes(dataBits, {'double','int8'}, ...
            {'real','integer','vector','binary'}, ...
            mfilename, 'Data input');

        % Columnize and expand to a [1 Nu] cell
        dataCell = repmat({int8(dataBits(:))}, 1, numUsers);
    else % NDP
        % Empty data check if not NDP
        coder.internal.errorIf(any(cfgFormat.PSDULength > 0) && isempty(dataBits), ...
            'wlan:wlanWaveformGenerator:NoData');

        dataCell = {int8(dataBits(:))};
    end
end

% Number of bits in a PSDU for a single packet (convert bytes to bits)
numPSDUBits = cfgFormat.PSDULength*8;

% Repeat to provide initial state(s) for all users and packets
scramInit = repmat(scramblerInit, 1, ...
    numUsers/size(scramblerInit, 2)); % For all users
pktScramInit = scramInit(mod((0:numPackets-1).', size(scramInit, 1))+1, :);

% Get the sampling rate of the waveform
if inDSSSMode % DSSS format
    sr = 11e6;
    numTxAnt = 1;
    giType = ''; % for codegen
    FFTLen = 0;  % for codegen
    info = dsssInfo(cfgFormat);
    numPktSamples = info.NumPPDUSamples;

    lstf = []; % for codegen
    lltf = []; % for codegen
    lsig = []; % for codegen
else % OFDM format
    chanBW = cfgFormat.ChannelBandwidth;
    sr = real(str2double(chanBW(4:end)))*1e6;
    if isa(cfgFormat, 'wlanNonHTConfig')
        giType = 'Long';  % Always
        FFTLen = 64;
    else    % For VHT/HT formats
        giType = cfgFormat.GuardInterval;
        num20 = real(str2double(chanBW(4:end)))/20;
        FFTLen = 64*num20;
    end
    if (strcmp(chanBW, 'CBW10') || strcmp(chanBW, 'CBW5'))
        numTxAnt = 1;  % override and set to 1 only, for 802.11j/p
    else
        numTxAnt = cfgFormat.NumTransmitAntennas;
    end
    numPktSamples = real(s.NumPPDUSamples); % real for codegen

    % Generate the legacy preamble fields for each PSDU
    lstf = wlanLSTF(cfgFormat);
    lltf = wlanLLTF(cfgFormat);
    lsig = wlanLSIG(cfgFormat);
end

if isa(cfgFormat, 'wlanVHTConfig') % VHT format
    % VHT Format
    vhtsiga = wlanVHTSIGA(cfgFormat);
    vhtstf = wlanVHTSTF(cfgFormat);
    vhtltf = wlanVHTLTF(cfgFormat);
    vhtsigb = wlanVHTSIGB(cfgFormat);
    preamble = [lstf; lltf; lsig; vhtsiga; vhtstf; vhtltf; vhtsigb];

elseif isa(cfgFormat, 'wlanHTConfig') % HT-MF format

    htSig = wlanHTSIG(cfgFormat);
    htstf = wlanHTSTF(cfgFormat);
    htltf = wlanHTLTF(cfgFormat);
    preamble = [lstf; lltf; lsig; htSig; htstf; htltf];

elseif isa(cfgFormat, 'wlanNonHTConfig')

    if strcmp(cfgFormat.Modulation, 'OFDM')
        preamble = [lstf; lltf; lsig];
    else % DSSS
        preamble = [wlanDSSSPreamble(cfgFormat); wlanDSSSHeader(cfgFormat)];
    end
end

if windowing
    % Calculate parameters for windowing

    % IdleSample offset due to windowing 
    wlength = 2*ceil(winTransitTime*sr/2);
    bLen = wlength/2; % Number of samples overlap at end of packet
    aLen = bLen-1;    % Number of samples overlap at start of packet
    
    windowedPktLength = numPktSamples+wlength-1;
else    
    % Define unused windowing variables for codegen
    wlength = 0;
    windowedPktLength = numPktSamples+wlength-1;
    aLen = 0;
    bLen = 0;
end

% Define a matrix of total simulation length
numIdleSamples = round(sr*idleTime);
pktWithIdleLength = numPktSamples+numIdleSamples;
txWaveform = complex(zeros(numPackets*pktWithIdleLength,numTxAnt));

for i = 1:numPackets
    % Extract PSDU for the current packet
    psdu = getPSDUForCurrentPacket(dataCell, numPSDUBits, i);

    % Generate the PSDU with the correct scrambler initial state
    if  isa(cfgFormat, 'wlanVHTConfig')
        if any(cfgFormat.APEPLength > 0)
            data = wlanVHTData(psdu, cfgFormat, pktScramInit(i, :));
        else  % NDP
            data = complex(zeros(0, cfgFormat.NumTransmitAntennas));
        end
    elseif isa(cfgFormat, 'wlanHTConfig') % HT-MF format
        if cfgFormat.PSDULength > 0                    
            data = wlanHTData(psdu{1}, cfgFormat, pktScramInit(i, :));
        else  % NDP or sounding packet
            data = complex(zeros(0, cfgFormat.NumTransmitAntennas));
        end
    elseif isa(cfgFormat, 'wlanNonHTConfig') % NonHT format
        if strcmp(cfgFormat.Modulation, 'OFDM')
            data = wlanNonHTData(psdu{1}, cfgFormat, pktScramInit(i, :));
        else % DSSS
            data = wlanDSSSData(psdu{1}, cfgFormat);
        end
    end

    % Construct packet from preamble and data
    packet = [preamble; data];                

    if windowing
        % Window each packet
         windowedPacket = wlanWindowing(packet, FFTLen, wlength, giType, ...
             size(preamble,1));

        % Overlap-add the windowed packets
        if numPackets==1 && numIdleSamples==0 % Only one packet which wraps     
            txWaveform = windowedPacket(aLen+(1:numPktSamples), :);
            % Overlap start of packet with end
            txWaveform(1:bLen, :) = txWaveform(1:bLen, :)+ ...
                windowedPacket(end-bLen+1:end, :);
            % Overlap end of packet with start
            txWaveform(end-aLen+1:end, :) = txWaveform(end-aLen+1:end, :)+ ...
                windowedPacket(1:aLen, :);
        else
            if i==1 % First packet (which wraps)
                % First packet wraps to end of waveform
                txWaveform(1:(numPktSamples+bLen), :) = ...
                    windowedPacket(1+aLen:end, :);
                txWaveform(end-aLen+1:end, :) = windowedPacket(1:aLen, :);
            elseif i==numPackets && numIdleSamples==0 % Last packet which wraps
                % Last packet wraps to start of waveform
                startIdx = (i-1)*pktWithIdleLength-aLen+1;
                txWaveform(startIdx:end, :) = txWaveform(startIdx:end, :)+ ...
                    windowedPacket(1:end-bLen, :);
                txWaveform(1:bLen,:) = txWaveform(1:bLen, :)+ ...
                    windowedPacket(end-bLen+1:end, :);
            else % Packet does not wrap
                % Normal windowing overlap between packets
                idx = (i-1)*pktWithIdleLength-aLen+(1:windowedPktLength);
                txWaveform(idx, :) = txWaveform(idx, :)+windowedPacket;
            end
       end
    else
        % Construct entire waveform
        txWaveform((i-1)*pktWithIdleLength+(1:numPktSamples), :) = packet;
    end
end
end

function psdu = getPSDUForCurrentPacket(dataCell, numPSDUBitsPerPacket, ...
    packetIdx)
    numUsers = length(dataCell); % == length(numPSDUBits)
    psdu = repmat({int8(1)}, 1, numUsers); % Cannot use cell(1, numUsers) for codegen
    for u = 1:numUsers
        idx = mod((packetIdx-1)*numPSDUBitsPerPacket(u)+...
            (0:numPSDUBitsPerPacket(u)-1).', length(dataCell{u})) + 1;
        psdu{u} = dataCell{u}(idx);
    end
end

% [EOF]
