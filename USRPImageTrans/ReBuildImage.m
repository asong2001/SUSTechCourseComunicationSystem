% 6. 图像重构 Reconstruct Image*
%     clc
%     clear
%     load('rxWaveform.mat')
%     rxWaveform2=[rxWaveform;rxWaveform];
     
    function ReBuildImage(rxBit,packetSeq)
    
    lengthMACheader = 256; % MPDU header length in bits
    rxBitMatrix = cell2mat(rxBit); 
    rxData = rxBitMatrix(lengthMACheader+1:end,1:numel(packetSeq)-1);
 
    startSeq = find(packetSeq==0);
    rxData = circshift(rxData,[0 -(startSeq(1)-1)]);% Order MAC fragments

%（1）计算误码率
%     % Perform bit error rate (BER) calculation
%     hBER = comm.ErrorRate;
%     err = step(hBER,double(rxData(:)),txData(1:length(reshape(rxData,[],1))));
%     fprintf('  \nBit Error Rate (BER):\n');
%     fprintf('          Bit Error Rate (BER) = %0.5f.\n',err(1));
%     fprintf('          Number of bit errors = %d.\n', err(2));
%     fprintf('    Number of transmitted bits = %d.\n\n',length(txData));

%（2）从接收的信号重构图像
    % Recreate image from received data
    fprintf('\nConstructing image from received data.\n');
    
    LtxImage=422280;  %--------------------------------------------->
    
    str = reshape(sprintf('%d',rxData(1:LtxImage)),8,[]).';
    decdata = uint8(bin2dec(str));
  
    imsize =[115   153     3];
    receivedImage = reshape(decdata,imsize);

%（3）打印图像
    % Plot received image
 %   if exist('imFig', 'var') && ishandle(imFig) % If Tx figure is open
 %       figure(imFig); subplot(212); 
 %   else
 %       figure; 
        subplot(212);
 %   end
        imshow(receivedImage);
    title(sprintf('Received Image'));


