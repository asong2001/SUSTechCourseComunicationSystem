% 接收过程分为3步：1. 捕获信号，2. 将波形解码成比特，3. 将比特还原成图像
clc
clear
%1. 捕获USRP信号
%(1)设置USRP参数；注意设置要参考发射机参数，特别是帧大小，即发射波形点数。
prmQPSKReceiver.USRPCenterFrequency = 900e6;
prmQPSKReceiver.USRPGain = 25;
prmQPSKReceiver.RxBufferedFrames=1;
prmQPSKReceiver.Fs = 5e6; % IQ 速率；
prmQPSKReceiver.USRPDecimationFactor = 100e6/prmQPSKReceiver.Fs; 
prmQPSKReceiver.FrameSize=202720;
prmQPSKReceiver.USRPFrameLength = ...
    prmQPSKReceiver.FrameSize*prmQPSKReceiver.RxBufferedFrames;

%(2) 构造接收机对象
radio = comm.SDRuReceiver(...
        'IPAddress',            '192.168.10.2', ...
        'CenterFrequency',      prmQPSKReceiver.USRPCenterFrequency, ...
        'Gain',                 prmQPSKReceiver.USRPGain, ...
        'DecimationFactor',     prmQPSKReceiver.USRPDecimationFactor, ...
        'FrameLength',          prmQPSKReceiver.USRPFrameLength, ...
        'OutputDataType',       'double');

%(3) 循环捕获，直到成功捕获数据包    
errorIndex=0;
while (true)
   % 从USRP读取IQ信号；
    [corruptSignal, len] = step(radio);

   % 如果未能成功读取我们希望的长度，报错
      if len < prmQPSKReceiver.USRPFrameLength 
         errorIndex = errorIndex+1;
         disp ( 'Not enough samples returned!' ) ;
         disp(errorIndex)
      else
          
         rxWaveform = corruptSignal;    
         break; % 否则，跳出循环，进一步做数据包解码操作。 
      end
end
release(radio);

%% 离线程序
load('rxWaveform.mat') %-------------------->这一行用于离线程序测试
rxWaveform2=[rxWaveform;rxWaveform];
[rxBit,packetSeq]=LoopPktDecode(rxWaveform2);
ReBuildImage(rxBit,packetSeq)
