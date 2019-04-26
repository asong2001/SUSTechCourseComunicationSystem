function [txImage,fData]=Image2Bits
% 1. ����
% �������Ĺ����ǽ�ͼ��ת��Ϊ������

% 2. ʵ�鲽�裺
% 2.1 ����Ƶ����ʾ����
% hsa = dsp.SpectrumAnalyzer( ...
%     'SpectrumType',    'Power density', ...
%     'SpectralAverages', 10, ...
%     'YLimits',         [-150 -50], ...
%     'Title',           'Received Baseband WLAN Signal Spectrum', ...
%     'YLabel',          'Power spectral density');
% 
% % 2.2 ��������ͼ����
% hcd = comm.ConstellationDiagram('Title','Equalized WLAN Symbols','ShowReferenceConstellation',false);
                            
%2.4 �������Ʋ���
% һ. ����ͼ�����ɶ���������
% ��. ���������źŴ����802.11a��ʽ���ź�
% ��. �����ϱ�Ƶ���������ݣ��������źŷ����ȥ

% һ. ����ͼ�����ɶ���������
%��1��ͼ���ļ���
fileTx = 'peppers.png';   % Image file name

%��2����ȡ�ļ���fData��һ����ά���飬��MATLAB�У�ͼ��ߴ緵��ֵ���һλ��ʾά��
fData = imread(fileTx);   % Read image data from file

%��������������������������������������������������������������������������������
%��һ���ֵ������ǽ�һ�����ͼ����С����
%��3������ߴ�任����
scale = 0.3;                % Image scaling factor

%��4��ԭʼͼ��ߴ磬����size(fData)������ֵ ans = 384   512     3
origSize = size(fData);   % Original input image size

%��5����Ҫ�����ͼ��ߴ� Calculate new image size
scaledSize = max(floor(scale.*origSize(1:2)),1); 

heightIx = min(round(((1:scaledSize(1))-0.5)./scale+0.5),origSize(1));
widthIx = min(round(((1:scaledSize(2))-0.5)./scale+0.5),origSize(2));

%��6�������������
fData = fData(heightIx,widthIx,:); % Resize image
%��������������������������������������������������������������������������������

%��7����ͼ��ĳߴ� size(fData) ans = 192   256     3
%imsize = size(fData);              % Store new image size

%��8��ת����8bit�޷��ŵĶ���������
binData = dec2bin(fData(:),8); 

%��9������������������
txImage = reshape((binData-'0').',1,[]).'; % Create binary stream

%��10����ʾ��Ҫ�����ͼ��
figure(1);
subplot(211); 
    imshow(fData);
    title('Transmitted Image');
subplot(212);
    title('Received image will appear here...');
    set(gca,'Visible','off');
    set(findall(gca, 'type', 'text'), 'visible', 'on');
hold on  
