function SimParams = sdru16QAMtransmitter_init(platform)
% Set simulation parameters
% SimParams = sdruqpsktransmitter_init

%   Copyright 2012-2015 The MathWorks, Inc.

switch platform
  case {'B200','B210'}
    SimParams.MasterClockRate = 20e6;  %Hz
    SimParams.Fs = 200e3; % Sample rate
  case {'X300','X310'}
    SimParams.MasterClockRate = 120e6; %Hz
    SimParams.Fs = 240e3; % Sample rate
  case {'N200/N210/USRP2'}
    SimParams.MasterClockRate = 100e6; %Hz
    SimParams.Fs = 200e3; % Sample rate
  otherwise
    error(message('sdru:examples:UnsupportedPlatform', ...
      platform))
end

% General simulation parameters
SimParams.Upsampling = 4; % Upsampling factor
SimParams.Ts = 1/SimParams.Fs; % Sample time
SimParams.FrameSize = 100; % Number of modulated symbols per frame

% Tx parameters
SimParams.BarkerLength = 13; % Number of Barker code symbols
SimParams.DataLength = (SimParams.FrameSize - SimParams.BarkerLength)*4; % Number of data payload bits per frame
SimParams.MessageLength = 112; % Number of message bits per frame, 7 ASCII characters
SimParams.FrameCount = 100;

SimParams.RxBufferedFrames = 10; % Received buffer length (in frames)
SimParams.RCFiltSpan = 10; % Filter span of Raised Cosine Tx Rx filters (in symbols)
SimParams.ScramblerBase = 2;
SimParams.ScramblerPolynomial = [1 1 1 0 1];
SimParams.ScramblerInitialConditions = [0 0 0 0];

% Generate square root raised cosine filter coefficients (required only for MATLAB example)
SimParams.SquareRootRaisedCosineFilterOrder = SimParams.Upsampling*SimParams.RCFiltSpan;
SimParams.RollOff = 0.5;

% Square root raised cosine transmit filter
hTxFilt = fdesign.interpolator(SimParams.Upsampling, ...
                'Square Root Raised Cosine', SimParams.Upsampling, ...
                'N,Beta', SimParams.SquareRootRaisedCosineFilterOrder, SimParams.RollOff);
hDTxFilt = design(hTxFilt, 'SystemObject', true);
SimParams.TransmitterFilterCoefficients = hDTxFilt.Numerator/2;

%SDRu transmitter parameters
SimParams.USRPCenterFrequency = 3.7e9;
SimParams.USRPGain = 25;
SimParams.USRPInterpolationFactor = SimParams.MasterClockRate/SimParams.Fs;
SimParams.USRPFrameLength = SimParams.Upsampling*SimParams.FrameSize*SimParams.RxBufferedFrames;

%Simulation Parameters
SimParams.FrameTime = SimParams.USRPFrameLength/SimParams.Fs;
SimParams.StopTime = 1000;
