function BER = run16QAMSystemUnderTest(prm16QAMTxRx, useScopes, printData)

% Copyright 2012-2014 The MathWorks, Inc.

%#codegen 

persistent hTx hChan hRx hScopes 
coder.extrinsic('create16QAMScopes','step16QAMScopes','release16QAMScopes')
if isempty(hTx)
    % Initialize the components
    % Create and configure the transmitter System object
    hTx = My16QAMTransmitter(...
        'UpsamplingFactor', prm16QAMTxRx.Upsampling, ...
        'MessageLength', prm16QAMTxRx.MessageLength, ...
        'TransmitterFilterCoefficients',prm16QAMTxRx.TransmitterFilterCoefficients, ...
        'DataLength', prm16QAMTxRx.DataLength, ...
        'ScramblerBase', prm16QAMTxRx.ScramblerBase, ...
        'ScramblerPolynomial', prm16QAMTxRx.ScramblerPolynomial, ...
        'ScramblerInitialConditions', prm16QAMTxRx.ScramblerInitialConditions);
    
    % Create and configure the AWGN channel System object
    hChan = My16QAMChannel('DelayType', prm16QAMTxRx.DelayType, ...
        'RaisedCosineFilterSpan', prm16QAMTxRx.RaisedCosineFilterSpan, ...
        'PhaseOffset', prm16QAMTxRx.PhaseOffset, ...
        'SignalPower', 1/prm16QAMTxRx.Upsampling, ...
        'FrameSize', prm16QAMTxRx.FrameSize, ...
        'UpsamplingFactor', prm16QAMTxRx.Upsampling, ...
        'EbNo', prm16QAMTxRx.EbNo, ...
        'BitsPerSymbol', prm16QAMTxRx.Upsampling/prm16QAMTxRx.Downsampling, ...
        'FrequencyOffset', prm16QAMTxRx.FrequencyOffset, ...
        'SampleRate', prm16QAMTxRx.Fs);

    % Create and configure the receiver System object
    hRx = My16QAMReceiver('DesiredAmplitude', 1/sqrt(prm16QAMTxRx.Upsampling), ...
        'ModulationOrder', prm16QAMTxRx.M, ...
        'DownsamplingFactor', prm16QAMTxRx.Downsampling, ...
        'CoarseCompFrequencyResolution', prm16QAMTxRx.CoarseCompFrequencyResolution, ...
        'PhaseRecoveryDampingFactor', prm16QAMTxRx.PhaseRecoveryDampingFactor, ...
        'PhaseRecoveryLoopBandwidth', prm16QAMTxRx.PhaseRecoveryLoopBandwidth, ...
        'TimingRecoveryDampingFactor', prm16QAMTxRx.TimingRecoveryDampingFactor, ...
        'TimingRecoveryLoopBandwidth', prm16QAMTxRx.TimingRecoveryLoopBandwidth, ...
        'TimingErrorDetectorGain', prm16QAMTxRx.TimingErrorDetectorGain, ...
        'PostFilterOversampling', prm16QAMTxRx.Upsampling/prm16QAMTxRx.Downsampling, ...
        'FrameSize', prm16QAMTxRx.FrameSize, ...
        'BarkerLength', prm16QAMTxRx.BarkerLength, ...
        'MessageLength', prm16QAMTxRx.MessageLength, ...
        'SampleRate', prm16QAMTxRx.Fs, ...
        'DataLength', prm16QAMTxRx.DataLength, ...
        'ReceiverFilterCoefficients', prm16QAMTxRx.ReceiverFilterCoefficients, ...
        'DescramblerBase', prm16QAMTxRx.ScramblerBase, ...
        'DescramblerPolynomial', prm16QAMTxRx.ScramblerPolynomial, ...
        'DescramblerInitialConditions', prm16QAMTxRx.ScramblerInitialConditions,...
        'PrintOption', printData);    
    
    if useScopes
        % Create the System object for plotting all the scopes
        hScopes = create16QAMScopes;
    end
end

% hRx.PrintOption = printData;

for count = 1:prm16QAMTxRx.FrameCount
    transmittedSignal = step(hTx); % Transmitter
    corruptSignal = step(hChan, transmittedSignal, count); % AWGN Channel
    [RCRxSignal,coarseCompBuffer, timingRecBuffer,BER] = step(hRx,corruptSignal); % Receiver
    if useScopes
        step16QAMScopes(hScopes,RCRxSignal,coarseCompBuffer, timingRecBuffer); % Plots all the scopes
    end
end
if isempty(coder.target)
    release(hTx);
    release(hChan);
    release(hRx);
end
if useScopes
     release16QAMScopes(hScopes);
end
