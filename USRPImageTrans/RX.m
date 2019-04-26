% ���չ��̷�Ϊ3����1. �����źţ�2. �����ν���ɱ��أ�3. �����ػ�ԭ��ͼ��
clc
clear
%1. ����USRP�ź�
%(1)����USRP������ע������Ҫ�ο�������������ر���֡��С�������䲨�ε�����
prmQPSKReceiver.USRPCenterFrequency = 900e6;
prmQPSKReceiver.USRPGain = 25;
prmQPSKReceiver.RxBufferedFrames=1;
prmQPSKReceiver.Fs = 5e6; % IQ ���ʣ�
prmQPSKReceiver.USRPDecimationFactor = 100e6/prmQPSKReceiver.Fs; 
prmQPSKReceiver.FrameSize=202720;
prmQPSKReceiver.USRPFrameLength = ...
    prmQPSKReceiver.FrameSize*prmQPSKReceiver.RxBufferedFrames;

%(2) ������ջ�����
radio = comm.SDRuReceiver(...
        'IPAddress',            '192.168.10.2', ...
        'CenterFrequency',      prmQPSKReceiver.USRPCenterFrequency, ...
        'Gain',                 prmQPSKReceiver.USRPGain, ...
        'DecimationFactor',     prmQPSKReceiver.USRPDecimationFactor, ...
        'FrameLength',          prmQPSKReceiver.USRPFrameLength, ...
        'OutputDataType',       'double');

%(3) ѭ������ֱ���ɹ��������ݰ�    
errorIndex=0;
while (true)
   % ��USRP��ȡIQ�źţ�
    [corruptSignal, len] = step(radio);

   % ���δ�ܳɹ���ȡ����ϣ���ĳ��ȣ�����
      if len < prmQPSKReceiver.USRPFrameLength 
         errorIndex = errorIndex+1;
         disp ( 'Not enough samples returned!' ) ;
         disp(errorIndex)
      else
          
         rxWaveform = corruptSignal;    
         break; % ��������ѭ������һ�������ݰ���������� 
      end
end
release(radio);

%% ���߳���
load('rxWaveform.mat') %-------------------->��һ���������߳������
rxWaveform2=[rxWaveform;rxWaveform];
[rxBit,packetSeq]=LoopPktDecode(rxWaveform2);
ReBuildImage(rxBit,packetSeq)
