function [PSDUData,numMSDUs]=Bits2PSDU(txImage)

%load('txImage.mat')

% 函数功能：创建MSDU->PSDU

% (1) MSDU字节数规定为4048个字节，不足的情况下补0.
   msduLength = 4048; % MSDU length in bytes

%（2）单个MSDU比特数
   msduBits = msduLength*8;

%（3）计算出需要分割的MSDU数目,本例子中一共有37个MSDU
   numMSDUs = ceil(length(txImage)/msduBits);

%（4）最后一个MSDU需要补0
   padZeros = msduBits-mod(length(txImage),msduBits);

%（5）需要发射数据
   txData = [txImage; zeros(padZeros,1)];

%（6）产生FCS
   generatorPolynomial = [32 26 23 22 16 12 11 10 8 7 5 4 2 1 0];
   fcsGen = comm.CRCGenerator(generatorPolynomial);
   fcsGen.InitialConditions = 1;
   fcsGen.DirectMethod = true;
   fcsGen.FinalXOR = 1;
 
%（7）将数据分块
   numFragment = 0;
%   bitsPerOctet = 8; 

%（8）MPDU头文件的比特数
   lengthMACheader = 256; % MPDU header length in bits

%（9）FCS帧校验比特数
   lengthFCS = 32;        % FCS length in bits

%（10）MPDU长度等于MAC头长度+MSDU比特+帧校验位,例如在本例子中，256+4048*8+32=1208864
   lengthMPDU = lengthMACheader+msduBits+lengthFCS; % MPDU length in bits

%（11）数据初始化
   PSDUData = zeros(lengthMPDU*numMSDUs,1);

%（12）形成MSDU数据包
for ind=0:numMSDUs-1

   %   构成MPDU
   % Extract bits for each MPDU
      frameBody = txData(ind*msduBits+1:msduBits*(ind+1),:);
   
   %   产生MPDU
   % Generate MPDU header bits
     mpduHeader = helperNonHTMACHeader(mod(numFragment, ...
                    16),mod(ind,4096));
   
   %   创建携带MAC头，帧体和FCS的MPDU数据
   % Create MPDU with header, body and FCS             
      psdu = step(fcsGen,[mpduHeader;frameBody]);
   
   %   传递到PSDU
   % Concatenate PSDUs for waveform generation
      PSDUData(lengthMPDU*ind+1:lengthMPDU*(ind+1)) = psdu;
   
end