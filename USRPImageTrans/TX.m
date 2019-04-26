% Transmitter 
% 分为两步：1. 802.11a波形产生；2. 利用USRP将波形发射出去；

clc
clear

%1. 802.11a波形产生，3小步来实现：
%（1）图像转换为比特流（信源编码）；
%（2）比特封装为PSDU，增加了MAC头文件；
%（3）将PSDU封装为PPDU，增加了PLCP头文件；

%(1)图像转换为比特流（信源编码）；
[txImage,fData]=Image2Bits;

%(2)比特封装为PSDU，增加了MAC头文件；
[PSDUData,numPSDUs]=Bits2PSDU(txImage);

%(3)将PSDU封装为PPDU，增加了PLCP头文件；
[txWaveform]=PSDU2PPDU(PSDUData,numPSDUs);
length(txWaveform)

%2. 利用USRP将波形发射出去；3小步来实现：
%(1)设置USRP参数；
prmQPSKTransmitter.USRPCenterFrequency = 900e6;
prmQPSKTransmitter.USRPGain = 25;
prmQPSKTransmitter.RxBufferedFrames=1;
prmQPSKTransmitter.Fs = 5e6; % IQ 速率；
prmQPSKTransmitter.USRPInterpolation = 100e6/prmQPSKTransmitter.Fs; % USRP上采样因子是100M/5M=20
prmQPSKTransmitter.FrameSize=length(txWaveform);
prmQPSKTransmitter.USRPFrameLength = prmQPSKTransmitter.FrameSize*prmQPSKTransmitter.RxBufferedFrames;

%Simulation Parameters
prmQPSKTransmitter.FrameTime = prmQPSKTransmitter.USRPFrameLength/prmQPSKTransmitter.Fs;
prmQPSKTransmitter.StopTime = 1000;

%（2）构造USRP发射机对象
    ThSDRu = comm.SDRuTransmitter('192.168.10.2', ...
        'CenterFrequency',        prmQPSKTransmitter.USRPCenterFrequency, ...
        'Gain',                   prmQPSKTransmitter.USRPGain, ...
        'InterpolationFactor',    prmQPSKTransmitter.USRPInterpolation);

%（3）循环发送
currentTime=0;
while currentTime < prmQPSKTransmitter.StopTime
    % Data transmission
    step(ThSDRu, txWaveform);
        
    % Update simulation time
    currentTime=currentTime+prmQPSKTransmitter.FrameTime
    
end

release(ThSDRu);
