
function BER = runQPSKSystemUnderTest(prmQPSKTxRx, useScopes, printData)
%

% Copyright 2012-2017 The MathWorks, Inc.

%#codegen

persistent qpskTx qpskChan qpskRx qpskScopes 
coder.extrinsic('createQPSKScopes','runQPSKScopes','releaseQPSKScopes')
if isempty(qpskTx)
    % Initialize the components
    % Create and configure the transmitter System object
    qpskTx = QPSKTransmitter(...
        'UpsamplingFactor',                     prmQPSKTxRx.Interpolation, ...
        'RolloffFactor',                        prmQPSKTxRx.RolloffFactor, ...
        'RaisedCosineFilterSpan',               prmQPSKTxRx.RaisedCosineFilterSpan, ...
        'MessageBits',                          prmQPSKTxRx.MessageBits, ...
        'MessageLength',                        prmQPSKTxRx.MessageLength, ...
        'NumberOfMessage',                      prmQPSKTxRx.NumberOfMessage, ...
        'ScramblerBase',                        prmQPSKTxRx.ScramblerBase, ...
        'ScramblerPolynomial',                  prmQPSKTxRx.ScramblerPolynomial, ...
        'ScramblerInitialConditions',           prmQPSKTxRx.ScramblerInitialConditions);
    
    % Create and configure the AWGN channel System object
    qpskChan = QPSKChannel(...
        'DelayType',                            prmQPSKTxRx.DelayType, ...
        'DelayStepSize',                        0.0125*prmQPSKTxRx.Interpolation, ...
        'DelayMaximum',                         2*prmQPSKTxRx.Interpolation, ...
        'DelayMinimum',                         0.1, ...
        'RaisedCosineFilterSpan',               prmQPSKTxRx.RaisedCosineFilterSpan, ...
        'PhaseOffset',                          prmQPSKTxRx.PhaseOffset, ...
        'SignalPower',                          1/prmQPSKTxRx.Interpolation, ...
        'InterpolationFactor',                  prmQPSKTxRx.Interpolation, ...
        'EbNo',                                 prmQPSKTxRx.EbNo, ...
        'BitsPerSymbol',                        log2(prmQPSKTxRx.ModulationOrder), ...
        'FrequencyOffset',                      prmQPSKTxRx.FrequencyOffset, ...
        'SampleRate',                           prmQPSKTxRx.Fs);

    % Create and configure the receiver System object
    qpskRx = QPSKReceiver(...
        'ModulationOrder',                      prmQPSKTxRx.ModulationOrder, ...
        'SampleRate',                           prmQPSKTxRx.Fs, ...
        'DecimationFactor',                     prmQPSKTxRx.Decimation, ...
        'FrameSize',                            prmQPSKTxRx.FrameSize, ...
        'HeaderLength',                         prmQPSKTxRx.HeaderLength, ...
        'NumberOfMessage',                      prmQPSKTxRx.NumberOfMessage, ...
        'PayloadLength',                        prmQPSKTxRx.PayloadLength, ...
        'DesiredPower',                         prmQPSKTxRx.DesiredPower, ...
        'AveragingLength',                      prmQPSKTxRx.AveragingLength, ...
        'MaxPowerGain',                         prmQPSKTxRx.MaxPowerGain, ...
        'RolloffFactor',                        prmQPSKTxRx.RolloffFactor, ...
        'RaisedCosineFilterSpan',               prmQPSKTxRx.RaisedCosineFilterSpan, ...
        'InputSamplesPerSymbol',                prmQPSKTxRx.Interpolation, ...
        'MaximumFrequencyOffset',               prmQPSKTxRx.MaximumFrequencyOffset, ...
        'PostFilterOversampling',               prmQPSKTxRx.Interpolation/prmQPSKTxRx.Decimation, ...
        'PhaseRecoveryLoopBandwidth',           prmQPSKTxRx.PhaseRecoveryLoopBandwidth, ...
        'PhaseRecoveryDampingFactor',           prmQPSKTxRx.PhaseRecoveryDampingFactor, ...
        'TimingRecoveryDampingFactor',          prmQPSKTxRx.TimingRecoveryDampingFactor, ...
        'TimingRecoveryLoopBandwidth',          prmQPSKTxRx.TimingRecoveryLoopBandwidth, ...
        'TimingErrorDetectorGain',              prmQPSKTxRx.TimingErrorDetectorGain, ...
        'PreambleDetectorThreshold',            prmQPSKTxRx.PreambleDetectorThreshold, ...    
        'DescramblerBase',                      prmQPSKTxRx.ScramblerBase, ...
        'DescramblerPolynomial',                prmQPSKTxRx.ScramblerPolynomial, ...
        'DescramblerInitialConditions',         prmQPSKTxRx.ScramblerInitialConditions,...
        'BerMask',                              prmQPSKTxRx.BerMask, ...
        'PrintOption',                          printData);

    if useScopes
        % Create the System object for plotting all the scopes
        sampleRate = prmQPSKTxRx.Rsym*prmQPSKTxRx.Interpolation/prmQPSKTxRx.Decimation;
        qpskScopes = createQPSKScopes(sampleRate);
    end
end

qpskRx.PrintOption = printData;

for count = 1:prmQPSKTxRx.TotalFrame
    transmittedSignal = qpskTx();                                           % Transmitter
    rcvdSignal = qpskChan(transmittedSignal, count);                        % AWGN Channel
    [RCRxSignal, timingRecSignal, freqRecSignal, BER] = qpskRx(rcvdSignal); % Receiver
    if useScopes
        runQPSKScopes(qpskScopes, rcvdSignal, RCRxSignal, timingRecSignal, freqRecSignal); % Plots all the scopes
    end
end
if isempty(coder.target)
    release(qpskTx);
    release(qpskChan);
    release(qpskRx);
end
if useScopes
     releaseQPSKScopes(qpskScopes);
end
