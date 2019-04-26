

 function [txWaveform]=PSDU2PPDU(PSDUData,numPSDUs)
%�������ܣ���PSDU���ݷ�װ��PPDU�������ɲ���

%load('PSDUdata.mat')
%��1������WLAN���ݰ�
    nonHTcfg = wlanNonHTConfig;         % Create packet configuration
    nonHTcfg.ChannelBandwidth='CBW5';
    
%��2�����Ʒ�ʽ��64QAM���ŵ���������Ϊ2/3
    nonHTcfg.MCS = 4;           % Modulation: 64QAM Rate: 2/3

%��3������������Ŀ
    nonHTcfg.NumTransmitAntennas = 1;   % Number of transmit antenna

%��4��ռ�ô���ʵ���У��ܵ�USRP�������ƣ������Ȳ���5MHz����
%     chanBW = nonHTcfg.ChannelBandwidth;

%��5������PSDU����(�ֽ�)
%     bitsPerOctet = 8;
%     nonHTcfg.PSDULength = lengthMPDU/bitsPerOctet; % Set the PSDU length
    nonHTcfg.PSDULength = 4084;
    
%��6����ʼ������
    scramblerInitialization = randi([1 127],numPSDUs,1);
 
%��7����������NonHT���ݰ�
    txWaveform = wlanWaveformGenerator(PSDUData,nonHTcfg, ...
     'NumPackets',numPSDUs,'IdleTime',80e-6, ...
     'ScramblerInitialization',scramblerInitialization);
 
% % %��8���ز������ο����Ȳ���
%     fs = helperSampleRate(chanBW); 
%     osf = 1.5;                     
%     txWaveform  = resample(txWaveform,fs*osf,fs);
    fprintf('\nGenerating WLAN transmit waveform:\n')

%��9����һ���ź�
    powerScaleFactor = 0.8;
    txWaveform = txWaveform.*(1/max(abs(txWaveform))*powerScaleFactor);
    
    figure
    subplot(211)
    plot(real(txWaveform(1:160)))
    subplot(212)
    plot(imag(txWaveform(1:160)))
    