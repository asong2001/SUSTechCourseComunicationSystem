clc;
clear;

%% USRP 参数
prmQPSKReceiver.USRPCenterFrequency = 900e6;
prmQPSKReceiver.USRPGain = 25;
prmQPSKReceiver.RxBufferedFrames=1;
prmQPSKReceiver.Fs = 5e6; % IQ 速率；
prmQPSKReceiver.USRPDecimationFactor = 100e6 / prmQPSKReceiver.Fs;
prmQPSKReceiver.FrameSize = 307200;
prmQPSKReceiver.USRPFrameLength = ...
    prmQPSKReceiver.FrameSize*prmQPSKReceiver.RxBufferedFrames;


%% receiver object
radio = comm.SDRuReceiver(...
    '192.168.10.2',...
    'CenterFrequency',  prmQPSKReceiver.USRPCenterFrequency,...
    'Gain',             prmQPSKReceiver.USRPGain,...
    'DecimationFactor', prmQPSKReceiver.USRPDecimationFactor,...
    'FrameLength',      prmQPSKReceiver.USRPFrameLength,...
    'OutputDataType' , 'double');


%% read data
errorIndex = 0;
while (true)
    [corruptSignal, len] = step(radio);

    if len < prmQPSKReceiver.USRPFrameLength
        errorIndex = errorIndex + 1;
        disp('Not enough samples returned');
        disp(errorIndex);
    else
        rxWaveform = corruptSignal;
        break;
    end
    
end

%% decode
disp('Constructing image from received data');
decodeRxDataStream = [decodeRxDataStream; decodeRxDataStream];
str = reshape(sprintf('%d', decodeRxDataStream(1:42280)),8,[]).';
decdata = uint8(bi2dex(str));
imsize = [115 153 3];
receuvedImage = reshape(decdata, imsize);
