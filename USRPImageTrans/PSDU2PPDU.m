

 function [txWaveform]=PSDU2PPDU(PSDUData,numPSDUs)
%函数功能：将PSDU数据封装成PPDU，并生成波形

%load('PSDUdata.mat')
%（1）创建WLAN数据包
    nonHTcfg = wlanNonHTConfig;         % Create packet configuration
    nonHTcfg.ChannelBandwidth='CBW5';
    
%（2）调制方式：64QAM，信道编码速率为2/3
    nonHTcfg.MCS = 4;           % Modulation: 64QAM Rate: 2/3

%（3）发射天线数目
    nonHTcfg.NumTransmitAntennas = 1;   % Number of transmit antenna

%（4）占用带宽，实验中，受到USRP速率限制，可以先采用5MHz带宽
%     chanBW = nonHTcfg.ChannelBandwidth;

%（5）设置PSDU长度(字节)
%     bitsPerOctet = 8;
%     nonHTcfg.PSDULength = lengthMPDU/bitsPerOctet; % Set the PSDU length
    nonHTcfg.PSDULength = 4084;
    
%（6）初始化扰码
    scramblerInitialization = randi([1 127],numPSDUs,1);
 
%（7）产生基带NonHT数据包
    txWaveform = wlanWaveformGenerator(PSDUData,nonHTcfg, ...
     'NumPackets',numPSDUs,'IdleTime',80e-6, ...
     'ScramblerInitialization',scramblerInitialization);
 
% % %（8）重采样波形可以先不用
%     fs = helperSampleRate(chanBW); 
%     osf = 1.5;                     
%     txWaveform  = resample(txWaveform,fs*osf,fs);
    fprintf('\nGenerating WLAN transmit waveform:\n')

%（9）归一化信号
    powerScaleFactor = 0.8;
    txWaveform = txWaveform.*(1/max(abs(txWaveform))*powerScaleFactor);
    
    figure
    subplot(211)
    plot(real(txWaveform(1:160)))
    subplot(212)
    plot(imag(txWaveform(1:160)))
    