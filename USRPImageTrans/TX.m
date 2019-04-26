% Transmitter 
% ��Ϊ������1. 802.11a���β�����2. ����USRP�����η����ȥ��

clc
clear

%1. 802.11a���β�����3С����ʵ�֣�
%��1��ͼ��ת��Ϊ����������Դ���룩��
%��2�����ط�װΪPSDU��������MACͷ�ļ���
%��3����PSDU��װΪPPDU��������PLCPͷ�ļ���

%(1)ͼ��ת��Ϊ����������Դ���룩��
[txImage,fData]=Image2Bits;

%(2)���ط�װΪPSDU��������MACͷ�ļ���
[PSDUData,numPSDUs]=Bits2PSDU(txImage);

%(3)��PSDU��װΪPPDU��������PLCPͷ�ļ���
[txWaveform]=PSDU2PPDU(PSDUData,numPSDUs);
length(txWaveform)

%2. ����USRP�����η����ȥ��3С����ʵ�֣�
%(1)����USRP������
prmQPSKTransmitter.USRPCenterFrequency = 900e6;
prmQPSKTransmitter.USRPGain = 25;
prmQPSKTransmitter.RxBufferedFrames=1;
prmQPSKTransmitter.Fs = 5e6; % IQ ���ʣ�
prmQPSKTransmitter.USRPInterpolation = 100e6/prmQPSKTransmitter.Fs; % USRP�ϲ���������100M/5M=20
prmQPSKTransmitter.FrameSize=length(txWaveform);
prmQPSKTransmitter.USRPFrameLength = prmQPSKTransmitter.FrameSize*prmQPSKTransmitter.RxBufferedFrames;

%Simulation Parameters
prmQPSKTransmitter.FrameTime = prmQPSKTransmitter.USRPFrameLength/prmQPSKTransmitter.Fs;
prmQPSKTransmitter.StopTime = 1000;

%��2������USRP���������
    ThSDRu = comm.SDRuTransmitter('192.168.10.2', ...
        'CenterFrequency',        prmQPSKTransmitter.USRPCenterFrequency, ...
        'Gain',                   prmQPSKTransmitter.USRPGain, ...
        'InterpolationFactor',    prmQPSKTransmitter.USRPInterpolation);

%��3��ѭ������
currentTime=0;
while currentTime < prmQPSKTransmitter.StopTime
    % Data transmission
    step(ThSDRu, txWaveform);
        
    % Update simulation time
    currentTime=currentTime+prmQPSKTransmitter.FrameTime
    
end

release(ThSDRu);
