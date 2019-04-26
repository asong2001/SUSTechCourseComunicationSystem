function [PSDUData,numMSDUs]=Bits2PSDU(txImage)

%load('txImage.mat')

% �������ܣ�����MSDU->PSDU

% (1) MSDU�ֽ����涨Ϊ4048���ֽڣ����������²�0.
   msduLength = 4048; % MSDU length in bytes

%��2������MSDU������
   msduBits = msduLength*8;

%��3���������Ҫ�ָ��MSDU��Ŀ,��������һ����37��MSDU
   numMSDUs = ceil(length(txImage)/msduBits);

%��4�����һ��MSDU��Ҫ��0
   padZeros = msduBits-mod(length(txImage),msduBits);

%��5����Ҫ��������
   txData = [txImage; zeros(padZeros,1)];

%��6������FCS
   generatorPolynomial = [32 26 23 22 16 12 11 10 8 7 5 4 2 1 0];
   fcsGen = comm.CRCGenerator(generatorPolynomial);
   fcsGen.InitialConditions = 1;
   fcsGen.DirectMethod = true;
   fcsGen.FinalXOR = 1;
 
%��7�������ݷֿ�
   numFragment = 0;
%   bitsPerOctet = 8; 

%��8��MPDUͷ�ļ��ı�����
   lengthMACheader = 256; % MPDU header length in bits

%��9��FCS֡У�������
   lengthFCS = 32;        % FCS length in bits

%��10��MPDU���ȵ���MACͷ����+MSDU����+֡У��λ,�����ڱ������У�256+4048*8+32=1208864
   lengthMPDU = lengthMACheader+msduBits+lengthFCS; % MPDU length in bits

%��11�����ݳ�ʼ��
   PSDUData = zeros(lengthMPDU*numMSDUs,1);

%��12���γ�MSDU���ݰ�
for ind=0:numMSDUs-1

   %   ����MPDU
   % Extract bits for each MPDU
      frameBody = txData(ind*msduBits+1:msduBits*(ind+1),:);
   
   %   ����MPDU
   % Generate MPDU header bits
     mpduHeader = helperNonHTMACHeader(mod(numFragment, ...
                    16),mod(ind,4096));
   
   %   ����Я��MACͷ��֡���FCS��MPDU����
   % Create MPDU with header, body and FCS             
      psdu = step(fcsGen,[mpduHeader;frameBody]);
   
   %   ���ݵ�PSDU
   % Concatenate PSDUs for waveform generation
      PSDUData(lengthMPDU*ind+1:lengthMPDU*(ind+1)) = psdu;
   
end