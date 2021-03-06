    %general simu params
SimParams.M=4;
SimParams.Upsampling=4;
SimParams.UpsamplingFactor = 4;
SimParams.Downsampling=2;
SimParams.DownsamplingFactor=2;
SimParams.Fs=2e5;
SimParams.Ts=1/SimParams.Fs;
SimParams.FrameSize=174;

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

%Tx params
SimParams.BarkerLength=26;
SimParams.DataLength=(SimParams.FrameSize-SimParams.BarkerLength)*2;
SimParams.ScramblerBase=2;
SimParams.ScramblerPolynomial=[1 1 1 0 1];
SimParams.ScramblerInitialConditions=[0 0 0 0];
% SimParams.sBit=sBit;
SimParams.RxBufferedFrames=10;
SimParams.RaisedCosineFilterSpan=10;
SimParams.MessageLength=112;
SimParams.FrameCount=100;

%channel params
SimParams.PhaseOffset=0;
SimParams.EbNo=2;
SimParams.FrequencyOffset=0;
SimParams.DelayType='Triangle';

%Rx params
SimParams.CoarseCompFrequencyResolution=25;

K=1;
A=1/sqrt(2);
SimParams.PhaseRecoveryLoopBandwidth=0.01;
SimParams.PhaseRecoveryDampingFactor=1;
SimParams.TimingRecoveryLoopBandwidth=0.01;
SimParams.TimingRecoveryDampingFactor=1;
SimParams.TimingErrorDetectorGain=2.7*2*K*A^2+2.7*2*K*A^2;

%QPSK mod Barker code header
BarkerCode=[1;1;1;1;1;-1;-1;1;1;-1;1;-1;1;1;1;1;1;1;-1;-1;1;1;-1;1;-1;1];
SimParams.ModulatedHeader=sqrt(2)/2*(-1-1i)*BarkerCode;

%generate square root raised cosine filter coefficients
SimParams.Rolloff=0.5;

%square root raised cosine transmit filter
SimParams.TransmitterFilterCoefficients=...
    rcosdesign(SimParams.Rolloff,SimParams.RaisedCosineFilterSpan,...
    SimParams.Upsampling);

%Squre root raised cosine receive filter
SimParams.ReceiverFilterCoefficients=...
    rcosdesign(SimParams.Rolloff,SimParams.RaisedCosineFilterSpan,...
    SimParams.Upsampling);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prmQPSKTxRx=SimParams;
printData=true;
useScopes=false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialize components
%create configure of Tx system object
hTx=QPSKTransmitterR(...
    'UpsamplingFactor',prmQPSKTxRx.Upsampling,...
    'MessageLength',prmQPSKTxRx.MessageLength,...
    'TransmitterFilterCoefficients',prmQPSKTxRx.TransmitterFilterCoefficients,...
    'DataLength',prmQPSKTxRx.DataLength,...
    'ScramblerBase',prmQPSKTxRx.ScramblerBase,...
    'ScramblerPolynomial',prmQPSKTxRx.ScramblerPolynomial,...
    'ScramblerInitialConditions',prmQPSKTxRx.ScramblerInitialConditions,...
    'LdpcNewH',prmQPSKTxRx.LdpcNewH,...
    'LdpcU',prmQPSKTxRx.LdpcU,...
    'LdpcL',prmQPSKTxRx.LdpcL);

%create configure of AWGN channel system object
hChan=QPSKChannelR('DelayType',prmQPSKTxRx.DelayType,...
    'RaisedCosineFilterSpan',prmQPSKTxRx.RaisedCosineFilterSpan,...
    'PhaseOffset',prmQPSKTxRx.PhaseOffset,...
    'SignalPower',1/prmQPSKTxRx.Upsampling,...
    'FrameSize',prmQPSKTxRx.FrameSize,...
    'UpsamplingFactor',prmQPSKTxRx.UpsamplingFactor,...
    'EbNo',prmQPSKTxRx.EbNo,...
    'BitsPerSymbol',prmQPSKTxRx.Upsampling/prmQPSKTxRx.Downsampling,...
    'FrequencyOffset',prmQPSKTxRx.FrequencyOffset,...
    'SampleRate',prmQPSKTxRx.Fs);

%create configure of Rx system object
hRx=QPSKReceiverR('DesiredAmplitude',1/sqrt(prmQPSKTxRx.Upsampling),...
    'ModulationOrder',prmQPSKTxRx.M,...
    'DownsamplingFactor',prmQPSKTxRx.DownsamplingFactor,...
    'CoarseFrequencyResolution',prmQPSKTxRx.CoarseCompFrequencyResolution,...
    'PhaseRecoveryDampingFactor',prmQPSKTxRx.PhaseRecoveryDampingFactor,...
    'PhaseRecoveryLoopBandwidth',prmQPSKTxRx.PhaseRecoveryLoopBandwidth,...
    'TimingRecoveryLoopBandwidth',prmQPSKTxRx.TimingRecoveryLoopBandwidth,...
    'TimingRecoveryDampingFactor',prmQPSKTxRx.TimingRecoveryDampingFactor,...
    'TimingErrorDetectorGain',prmQPSKTxRx.TimingErrorDetectorGain,...
    'PostFilterOversampling',prmQPSKTxRx.Upsampling/prmQPSKTxRx.Downsampling,...
    'FrameSize',prmQPSKTxRx.FrameSize,...
    'BarkerLength',prmQPSKTxRx.BarkerLength,...
    'MessageLength',prmQPSKTxRx.MessageLength,...
    'SampleRate',prmQPSKTxRx.Fs,...
    'DataLength',prmQPSKTxRx.DataLength,...
    'ReceiverFilterCoefficients',prmQPSKTxRx.ReceiverFilterCoefficients,...
    'DescramblerBase',prmQPSKTxRx.ScramblerBase,...
    'DescramblerPolynomial',prmQPSKTxRx.ScramblerPolynomial,...
    'DescramblerInitialConditions',prmQPSKTxRx.ScramblerInitialConditions,...
    'LdpcNewH',prmQPSKTxRx.LdpcNewH,...
    'LdpcIteration',prmQPSKTxRx.LdpcIteration,...
    'LdpcU',prmQPSKTxRx.LdpcU,...
    'LdpcL',prmQPSKTxRx.LdpcL,...
    'PrintOption',printData);

if useScopes
    hScopes=createQPSKScopes;
end

hRx.PrintOption=printData;

for count=1:prmQPSKTxRx.FrameCount
    [transmittedSignal]=step(hTx);
    cooruptSignal=step(hChan,transmittedSignal,0);
    [RCRxSignal,coarseCompBuffer,timingRecBuffer,BER]=step(hRx,cooruptSignal);
%     if useScopes
%         step16QAMScopes(hScopes,RCRxSignal,coarseCompBuffer, timingRecBuffer); % Plots all the scopes
%     end
    if useScopes
        runQPSKScopes(hScopes,RCRxSignal,coarseCompBuffer, timingRecBuffer); % Plots all the scopes
    end

end

if useScopes
    release16QAMScopes(hScopes);
end

fprintf('Error rate = %f.\n',BER(1));
fprintf('Number of detected errors = %d.\n',BER(2));
fprintf('Total number of compared samples = %d.\n',BER(3));

