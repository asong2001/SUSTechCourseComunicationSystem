clc;
clear;

% LTE Transmitter

%% Read image data and convert
% *Prepare Image File 图像->比特流
% ----------图片转化部分------------ % 
% fileTx = 'peppers.png';            % Image file name
fileTx = 'robot.png';            % Image file name
fData = imread(fileTx);            % Read image data from file
fData = imresize(fData, [115 153]);
scale = 0.3;                       % Image scaling factor
origSize = size(fData);            % Original input image size
scaledSize = max(floor(scale.*origSize(1:2)),1); % Calculate new image size
heightIx = min(round(((1:scaledSize(1))-0.5)./scale+0.5),origSize(1));
widthIx = min(round(((1:scaledSize(2))-0.5)./scale+0.5),origSize(2));
fData = fData(heightIx,widthIx,:); % Resize image
imsize = size(fData);              % Store new image size
binData = dec2bin(fData(:),8);     % Convert to 8 bit unsigned binary
trData = reshape((binData-'0').',1,[]).'; % Create binary stream

%% LTE frame data
[eNodeBOutput, txGrid, rmc] = lteRmCDLTool(rmc, trData);

txWavaform = eNodeBOutput;
powerScaleFactor = 0.8;
txWavaform = txWavaform.*(1/max(abs(txWavaform)) * powerScaleFactor);

disp(['Wavaform size', size(txWavaform)]);

%% USRP 参数
prmQPSKTransmitter.USRPCenterFrequency = 900e6;
prmQPSKTransmitter.USRPGain = 25;
prmQPSKTransmitter.RxBufferedFrames=1;
prmQPSKTransmitter.Fs = 5e6; % IQ 速率；
prmQPSKTransmitter.USRPInterplation = 100e6 / prmQPSKTransmitter.Fs;
prmQPSKTransmitter.FrameSize = length(txWavaform);
prmQPSKTransmitter.USRPFrameLength = ...
    prmQPSKTransmitter.FrameSize*prmQPSKTransmitter.RxBufferedFrames;

%% simulation parameters
prmQPSKTransmitter.FrameTime = ...
    prmQPSKTransmitter.USRPFrameLength / prmQPSKTransmitter.Fs;
prmQPSKTransmitter.StopTime = 1000;

%% object of USRP
ThSDRu = comm.SDRuTransmitter('192.168.10.2',...
    'CenterFrequency', prmQPSKTransmitter.USRPCenterFrequency,...
    'Gain',             prmQPSKTransmitter.USRPGain,...
    'InterpolationFactor', prmQPSKTransmitter.USRPInterplation);

%% USRP transmitting
currentTime = 0;
while currentTime < prmQPSKTransmitter.StopTime
    step(ThSDRu, txWavaform);

    currentTime = currentTime + prmQPSKTransmitter.FrameTime
    
end