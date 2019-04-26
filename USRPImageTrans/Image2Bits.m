function [txImage,fData]=Image2Bits
% 1. 概述
% 本函数的功能是将图像转换为比特流

% 2. 实验步骤：
% 2.1 配置频谱显示工具
% hsa = dsp.SpectrumAnalyzer( ...
%     'SpectrumType',    'Power density', ...
%     'SpectralAverages', 10, ...
%     'YLimits',         [-150 -50], ...
%     'Title',           'Received Baseband WLAN Signal Spectrum', ...
%     'YLabel',          'Power spectral density');
% 
% % 2.2 配置星座图工具
% hcd = comm.ConstellationDiagram('Title','Equalized WLAN Symbols','ShowReferenceConstellation',false);
                            
%2.4 发射机设计步骤
% 一. 导入图像，生成二进制码流
% 二. 将二进制信号打包成802.11a格式的信号
% 三. 利用上变频和连续传递，将基带信号发射出去

% 一. 导入图像，生成二进制码流
%（1）图像文件名
fileTx = 'peppers.png';   % Image file name

%（2）读取文件，fData是一个三维数组，在MATLAB中，图像尺寸返回值最后一位表示维度
fData = imread(fileTx);   % Read image data from file

%————————————————————————————————————————
%这一部分的作用是将一幅大的图像缩小传输
%（3）定义尺寸变换因子
scale = 0.3;                % Image scaling factor

%（4）原始图像尺寸，例如size(fData)，返回值 ans = 384   512     3
origSize = size(fData);   % Original input image size

%（5）需要传输的图像尺寸 Calculate new image size
scaledSize = max(floor(scale.*origSize(1:2)),1); 

heightIx = min(round(((1:scaledSize(1))-0.5)./scale+0.5),origSize(1));
widthIx = min(round(((1:scaledSize(2))-0.5)./scale+0.5),origSize(2));

%（6）重新组合数据
fData = fData(heightIx,widthIx,:); % Resize image
%————————————————————————————————————————

%（7）新图像的尺寸 size(fData) ans = 192   256     3
%imsize = size(fData);              % Store new image size

%（8）转换成8bit无符号的二进制数据
binData = dec2bin(fData(:),8); 

%（9）创建二进制数据流
txImage = reshape((binData-'0').',1,[]).'; % Create binary stream

%（10）显示需要传输的图像
figure(1);
subplot(211); 
    imshow(fData);
    title('Transmitted Image');
subplot(212);
    title('Received image will appear here...');
    set(gca,'Visible','off');
    set(findall(gca, 'type', 'text'), 'visible', 'on');
hold on  
