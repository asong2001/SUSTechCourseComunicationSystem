clc;
clear all;

SimParams.MasterClockRate = 100e6; %Hz
SimParams.Fs = 2e5; % Sample rate
% General simulation parameters
SimParams.Upsampling = 4; % Upsampling factor
SimParams.Ts = 1/SimParams.Fs; % Sample time
SimParams.FrameSize = 174; % Number of modulated symbols per frame

% Tx parameters
SimParams.BarkerLength = 26; % Number of Barker code symbols
SimParams.DataLength = (SimParams.FrameSize - SimParams.BarkerLength)*2; % Number of data payload bits per frame
SimParams.MessageLength = 112; % Number of message bits per frame, 7 ASCII characters
SimParams.FrameCount = 100;

%channel coding params
SimParams.LdpcRow=148;
SimParams.LdpcColumn=296;
SimParams.Ldpcmethod=1;
SimParams.LdpcnoCycle=1;
SimParams.LdpcOnePerCol=3;
SimParams.LdpcStrategy=1;
SimParams.LdpcIteration=1;
SimParams.LdpcH=makeLdpc(SimParams.LdpcRow,...
                        SimParams.LdpcColumn,...
                        SimParams.Ldpcmethod,...
                        SimParams.LdpcnoCycle,...
                        SimParams.LdpcOnePerCol);
[SimParams.LdpcNewH,SimParams.LdpcU,SimParams.LdpcL]=makeParityChk(SimParams.LdpcH,SimParams.LdpcStrategy);


SimParams.RxBufferedFrames = 10; % Received buffer length (in frames)
SimParams.RaiseCosineGroupDelay = 10; % Group Delay of Raised Cosine Tx Rx filters (in symbols)
SimParams.ScramblerBase = 2;
SimParams.ScramblerPolynomial = [1 1 1 0 1];
SimParams.ScramblerInitialConditions = [0 0 0 0];

% Generate square root raised cosine filter coefficients (required only for MATLAB example)
SimParams.SquareRootRaisedCosineFilterOrder = 2*SimParams.Upsampling*SimParams.RaiseCosineGroupDelay;
SimParams.RollOff = 0.5;

% Square root raised cosine transmit filter
ThTxFilt = fdesign.interpolator(SimParams.Upsampling, ...
                'Square Root Raised Cosine', SimParams.Upsampling, ...
                'N,Beta', SimParams.SquareRootRaisedCosineFilterOrder, SimParams.RollOff);
ThDTxFilt = design(ThTxFilt);
SimParams.TransmitterFilterCoefficients = ThDTxFilt.Numerator/2;

%SDRu transmitter parameters
SimParams.USRPCenterFrequency = 4e9;
SimParams.USRPGain = 25;
SimParams.USRPInterpolationFactor = SimParams.MasterClockRate/SimParams.Fs;
SimParams.USRPFrameLength = SimParams.Upsampling*SimParams.FrameSize*SimParams.RxBufferedFrames;

%Simulation Parameters
SimParams.FrameTime = SimParams.USRPFrameLength/SimParams.Fs;
SimParams.StopTime = 1000000;

%%
prmQPSKTransmitter = SimParams

% Initialize the components
% Create and configure the transmitter System object
hTx = QPSKTransmitterR(...
    'UpsamplingFactor', prmQPSKTransmitter.Upsampling, ...
    'MessageLength', prmQPSKTransmitter.MessageLength, ...
    'TransmitterFilterCoefficients',prmQPSKTransmitter.TransmitterFilterCoefficients, ...
    'DataLength', prmQPSKTransmitter.DataLength, ...
    'ScramblerBase', prmQPSKTransmitter.ScramblerBase, ...
    'ScramblerPolynomial', prmQPSKTransmitter.ScramblerPolynomial, ...
    'ScramblerInitialConditions', prmQPSKTransmitter.ScramblerInitialConditions,...
    'LdpcNewH',prmQPSKTxRx.LdpcNewH,...
    'LdpcU',prmQPSKTxRx.LdpcU,...
    'LdpcL',prmQPSKTxRx.LdpcL);

ThSDRu = comm.SDRuTransmitter('192.168.10.2',...
    'CenterFrequency', prmQPSKTransmitter.USRPCenterFrequency,...
    'Gain', prmQPSKTransmitter.USRPGain,...
    'InterpolationFactor', prmQPSKTransmitter.USRPInterpolationFactor);

currentTime = 0;

%Transmission Process
while currentTime < prmQPSKTransmitter.StopTime
    % Bit generation, modulation and transmission filtering
    data = step(hTx);
    % Data transmission
    step(ThSDRu, data);
    % Update simulation time
    currentTime=currentTime+prmQPSKTransmitter.FrameTime
end

release(hTx);
release(ThSDRu);