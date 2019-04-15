function SimParams = commqpsktxrx_init
% Set simulation parameters

% Copyright 2011-2018 The MathWorks, Inc.

%% General simulation parameters
SimParams.ModulationOrder = 4;      % QPSK alphabet size
SimParams.Interpolation = 2;        % Interpolation factor
SimParams.Decimation = 1;           % Decimation factor
SimParams.Rsym = 5e4;               % Symbol rate in Hertz
SimParams.Tsym = 1/SimParams.Rsym;  % Symbol time in sec
SimParams.Fs   = SimParams.Rsym * SimParams.Interpolation; % Sample rate
SimParams.TotalFrame = 10;        % Simulate 1000 frames in total

%% Frame Specifications
% [BarkerCode*2 | 'Hello world 000\n' | 'Hello world 001\n' ...];
SimParams.BarkerCode      = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];     % Bipolar Barker Code
SimParams.BarkerLength    = length(SimParams.BarkerCode);
SimParams.HeaderLength    = SimParams.BarkerLength * 2;                   % Duplicate 2 Barker codes to be as a header
SimParams.Message         = 'Hello world';
SimParams.MessageLength   = strlength(SimParams.Message) + 5;                % 'Hello world 000\n'...
SimParams.NumberOfMessage = 20;                                           % Number of messages in a frame
SimParams.PayloadLength   = SimParams.NumberOfMessage * SimParams.MessageLength * 7; % 7 bits per characters
SimParams.FrameSize       = (SimParams.HeaderLength + SimParams.PayloadLength) ...
    / log2(SimParams.ModulationOrder);                                    % Frame size in symbols
SimParams.FrameTime       = SimParams.Tsym*SimParams.FrameSize;

%% Tx parameters
SimParams.RolloffFactor     = 0.5;                                          % Rolloff Factor of Raised Cosine Filter
SimParams.ScramblerBase     = 2;
SimParams.ScramblerPolynomial           = [1 1 1 0 1];
SimParams.ScramblerInitialConditions    = [0 0 0 0];
SimParams.RaisedCosineFilterSpan = 10; % Filter span of Raised Cosine Tx Rx filters (in symbols)

%% Channel parameters
SimParams.PhaseOffset       = 47;   % in degrees
SimParams.EbNo              = 0;   % in dB
SimParams.FrequencyOffset   = 5000; % Frequency offset introduced by channel impairments in Hertz
SimParams.DelayType         = 'Triangle'; % select the type of delay for channel distortion

%% MIMO Channel
SimParams.numTx = 2;         % Number of transmit antennas
SimParams.numRx = 2;         % Number of receive antennas
SimParams.Rs = 1e6;          % Sampling rate (Hz)
SimParams.tau = [0 2e-5];    % Path delays (sec)
SimParams.pdb = [0 -10];     % Average path gains (dB)
SimParams.maxDopp = 30;      % Maximum Doppler shift (Hz)
SimParams.numBits = 12000;   % Number of bits
SimParams.SNR = 6;           % Signal-to-noise ratio (dB)

%% Rx parameters
SimParams.DesiredPower                  = 2;            % AGC desired output power (in watts)
SimParams.AveragingLength               = 50;           % AGC averaging length
SimParams.MaxPowerGain                  = 20;           % AGC maximum output power gain
SimParams.MaximumFrequencyOffset        = 6e3;
% Look into model for details for details of PLL parameter choice. Refer equation 7.30 of "Digital Communications - A Discrete-Time Approach" by Michael Rice.
K = 1;
A = 1/sqrt(2);
SimParams.PhaseRecoveryLoopBandwidth    = 0.01;         % Normalized loop bandwidth for fine frequency compensation
SimParams.PhaseRecoveryDampingFactor    = 1;            % Damping Factor for fine frequency compensation
SimParams.TimingRecoveryLoopBandwidth   = 0.01;         % Normalized loop bandwidth for timing recovery
SimParams.TimingRecoveryDampingFactor   = 1;            % Damping Factor for timing recovery
% K_p for Timing Recovery PLL, determined by 2KA^2*2.7 (for binary PAM), 
% QPSK could be treated as two individual binary PAM, 
% 2.7 is for raised cosine filter with roll-off factor 0.5
SimParams.TimingErrorDetectorGain       = 2.7*2*K*A^2+2.7*2*K*A^2; 
SimParams.PreambleDetectorThreshold     = 20;

%% Message generation and BER calculation parameters
msgSet = zeros(100 * SimParams.MessageLength, 1); 
for msgCnt = 0 : 99
    msgSet(msgCnt * SimParams.MessageLength + (1 : SimParams.MessageLength)) = ...
        sprintf('%s %03d\n', SimParams.Message, msgCnt);
end
integerToBit = comm.IntegerToBit(7, 'OutputDataType', 'double');
SimParams.MessageBits = integerToBit(msgSet);

% For BER calculation masks
SimParams.BerMask = zeros(SimParams.NumberOfMessage * length(SimParams.Message) * 7, 1);
for i = 1 : SimParams.NumberOfMessage
    SimParams.BerMask( (i-1) * length(SimParams.Message) * 7 + ( 1: length(SimParams.Message) * 7) ) = ...
        (i-1) * SimParams.MessageLength * 7 + (1: length(SimParams.Message) * 7);
end
