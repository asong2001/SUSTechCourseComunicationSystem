function runSDRu16QAMTransmitter(prm16QAMTransmitter)
%#codegen

%   Copyright 2012-2014 The MathWorks, Inc.

persistent hTx radio
if isempty(hTx)
  % Initialize the components
  % Create and configure the transmitter System object
  hTx = QAMTransmitter(...
    'UpsamplingFactor', prm16QAMTransmitter.Upsampling, ...
    'MessageLength', prm16QAMTransmitter.MessageLength, ...
    'TransmitterFilterCoefficients',prm16QAMTransmitter.TransmitterFilterCoefficients, ...
    'DataLength', prm16QAMTransmitter.DataLength, ...
    'ScramblerBase', prm16QAMTransmitter.ScramblerBase, ...
    'ScramblerPolynomial', prm16QAMTransmitter.ScramblerPolynomial, ...
    'ScramblerInitialConditions', prm16QAMTransmitter.ScramblerInitialConditions);
  
  % Create and configure the SDRu System object. Set the SerialNum for B2xx
  % radios and IPAddress for X3xx, N2xx, and USRP2 radios. MasterClockRate
  % is not configurable for N2xx and USRP2 radios.
  switch prm16QAMTransmitter.Platform
    case {'B200','B210'}
      radio = comm.SDRuTransmitter(...
        'Platform',             prm16QAMTransmitter.Platform, ...
        'SerialNum',            prm16QAMTransmitter.Address, ...
        'MasterClockRate',      prm16QAMTransmitter.MasterClockRate, ...
        'CenterFrequency',      prm16QAMTransmitter.USRPCenterFrequency, ...
        'Gain',                 prm16QAMTransmitter.USRPGain, ...
        'InterpolationFactor',  prm16QAMTransmitter.USRPInterpolationFactor);
    case {'X300','X310'}
      radio = comm.SDRuTransmitter(...
        'Platform',             prm16QAMTransmitter.Platform, ...
        'IPAddress',            prm16QAMTransmitter.Address, ...
        'MasterClockRate',      prm16QAMTransmitter.MasterClockRate, ...
        'CenterFrequency',      prm16QAMTransmitter.USRPCenterFrequency, ...
        'Gain',                 prm16QAMTransmitter.USRPGain, ...
        'InterpolationFactor',  prm16QAMTransmitter.USRPInterpolationFactor);
    case {'N200/N210/USRP2'}
      radio = comm.SDRuTransmitter(...
        'Platform',             prm16QAMTransmitter.Platform, ...
        'IPAddress',            prm16QAMTransmitter.Address, ...
        'CenterFrequency',      prm16QAMTransmitter.USRPCenterFrequency, ...
        'Gain',                 prm16QAMTransmitter.USRPGain, ...
        'InterpolationFactor',  prm16QAMTransmitter.USRPInterpolationFactor);
  end
  
  currentTime = 0;
  
  %Transmission Process
  while currentTime < prm16QAMTransmitter.StopTime
    % Bit generation, modulation and transmission filtering
    data = step(hTx);
    % Data transmission
    step(radio, data);
    % Update simulation time
    currentTime=currentTime+prm16QAMTransmitter.FrameTime;
  end
  
  release(hTx);
  release(radio);
end