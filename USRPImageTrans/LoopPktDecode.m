   %clc
   %clear
   function [rxBit,packetSeq]=LoopPktDecode(rxWaveform)
   
% ��1�������ز������������в���Ҫ
%   rxWaveform = resample(burstCaptures,fs,fs*osf);
    %load('rxWaveform.mat');
    rxWaveformLen = size(rxWaveform,1);
    searchOffset = 0; % Offset from start of the waveform in samples
  
% ��2����С��������10��OFDM����
    nonHTcfg = wlanNonHTConfig;
    chanBW='CBW5';     %--------------------------------->����ط�Ҫ�޸�
    indLSTF = wlanFieldIndices(nonHTcfg,'L-STF'); 
    indLLTF = wlanFieldIndices(nonHTcfg,'L-LTF'); 
    indLSIG = wlanFieldIndices(nonHTcfg,'L-SIG');
    lstfLen = double(indLSTF(2)); % Number of samples in L-STF
    minPktLen = lstfLen*5;
    pktInd = 1;
    sr = helperSampleRate(chanBW); % Sampling rate
    %offsetLLTF = [];
    packetSeq = [];
    displayFlag = 1; % Flag to display the decoded information

%��3��ΪMPDU����FCS
    generatorPolynomial = [32 26 23 22 16 12 11 10 8 7 5 4 2 1 0]; 
    fcsDet = comm.CRCDetector(generatorPolynomial);
    fcsDet.InitialConditions = 1;
    fcsDet.DirectMethod = true;
    fcsDet.FinalXOR = 1;

%��4������EVM
    hEVM = comm.EVM('AveragingDimensions',[1 2 3]);
    hEVM.MaximumEVMOutputPort = true;


%��5������ѭ������
while (searchOffset + minPktLen) <= rxWaveformLen    
    
% ���ݰ���� Packet detect
    pktOffset = helperPacketDetect(rxWaveform(1+searchOffset:end,:),chanBW,0.8)-1;
 
% �������ݰ�ƫ�� Adjust packet offset
    pktOffset = searchOffset+pktOffset;
    if isempty(pktOffset) || (pktOffset+indLSIG(2)>rxWaveformLen)
        if pktInd==1
            disp('** No packet detected **');
        end
        break;
    end
 
% ��ȡNon-HT�򣬴�Ƶƫ����
    nonHT = rxWaveform(pktOffset+(indLSTF(1):indLSIG(2)),:);
    coarseFreqOffset = wlanCoarseCFOEstimate(nonHT,chanBW); 
    nonHT = helperFrequencyOffset(nonHT,sr,-coarseFreqOffset);

% ����ͬ�� Symbol timing synchronization
    offsetLLTF = helperSymbolTiming(nonHT,chanBW);
    
    if isempty(offsetLLTF)
        searchOffset = pktOffset+lstfLen;
        continue;
    end
    % Adjust packet offset
    pktOffset = pktOffset+offsetLLTF-double(indLLTF(1));

% ��ʱͬ�����
    if (pktOffset<0) || ((pktOffset+minPktLen)>rxWaveformLen) 
        searchOffset = pktOffset+lstfLen; 
        continue; 
    end
    fprintf('\nPacket-%d detected at index %d\n',pktInd,pktOffset+1);
  
    nonHT = rxWaveform(pktOffset+(indLSTF(1):indLSIG(2)),:);
    nonHT = helperFrequencyOffset(nonHT,sr,-coarseFreqOffset);

% ��Ƶƫ����
    lltf = nonHT(indLLTF(1):indLLTF(2),:);           % Extract L-LTF
    fineFreqOffset = wlanFineCFOEstimate(lltf,chanBW);
    nonHT = helperFrequencyOffset(nonHT,sr,-fineFreqOffset);
    cfoCorrection = coarseFreqOffset+fineFreqOffset; % Total CFO
 
% ����L-LTF���ŵ�����
    lltf = nonHT(indLLTF(1):indLLTF(2),:);
    demodLLTF = wlanLLTFDemodulate(lltf,chanBW);
    chanEstLLTF = wlanLLTFChannelEstimate(demodLLTF,chanBW);

% ��������
    noiseVarNonHT = helperNoiseEstimate(demodLLTF);

% �ָ���L-SIG��
    [recLSIGBits,failCheck] = wlanLSIGRecover( ...
           nonHT(indLSIG(1):indLSIG(2),:), ...
           chanEstLLTF, noiseVarNonHT,chanBW);
 
    if failCheck	
        fprintf('  L-SIG check fail \n');
        searchOffset = pktOffset+lstfLen; 
        continue; 
    else
        fprintf('  L-SIG check pass \n');
    end
 
% �ָ����ݰ�����
    [lsigMCS,lsigLen,rxSamples] = helperInterpretLSIG(recLSIGBits,sr);
 
    if (rxSamples+pktOffset)>length(rxWaveform)
        disp('** Not enough samples to decode packet **');
        break;
    end

% Ӧ��CFO�����������ݰ�
    % Apply CFO correction to the entire packet
    rxWaveform(pktOffset+(1:rxSamples),:) = ...
        helperFrequencyOffset(rxWaveform(pktOffset+(1:rxSamples),:),sr,-cfoCorrection);

% ��������Non-HT����
    % Create a receive Non-HT config object
    rxNonHTcfg = wlanNonHTConfig;
    rxNonHTcfg.MCS = lsigMCS;
    rxNonHTcfg.PSDULength = lsigLen;

% ��ȡ������ָʾ
    % Get the data field indices within a PPDU 
    indNonHTData = wlanFieldIndices(rxNonHTcfg,'NonHT-Data');


% �����ŵ����ƽ���ָ�PSDU����
    [rxPSDU,eqSym] = wlanNonHTDataRecover(rxWaveform(pktOffset+...
           (indNonHTData(1):indNonHTData(2)),:), ...
           chanEstLLTF,noiseVarNonHT,rxNonHTcfg);
 
% ��ʾ��ǰ����ͼ
  %   hcd = comm.ConstellationDiagram('Title','Equalized WLAN Symbols','ShowReferenceConstellation',false);
  %   step(hcd,reshape(eqSym,[],1)); % Current constellation 
  %   release(hcd); % Release previous constellation plot
 
     refSym = helperClosestConstellationPoint(eqSym,rxNonHTcfg);
    [evm.RMS,evm.Peak] = step(hEVM,refSym,eqSym);

% ��MACͷ���Ƴ�FCS
    % Remove FCS from MAC header and frame body
    [rxBit{pktInd},crcCheck] = step(fcsDet,double(rxPSDU)); 
 
    if ~crcCheck
         disp('  MAC CRC check pass');
    else
         disp('  MAC CRC check fail');
    end
 
% ����MAC��Ϣ
    [mac,packetSeq(pktInd)] = helperNonHTMACHeaderDecode(rxBit{pktInd}); 

% ��ʾ������
    % Display decoded information
    if displayFlag
        fprintf('  Estimated CFO: %5.1f Hz\n\n',cfoCorrection); 
 
        disp('  Decoded L-SIG contents: ');
        fprintf(' MCS: %d\n',lsigMCS);
        fprintf(' Length: %d\n',lsigLen);
        fprintf(' Number of samples in packet: %d\n\n',rxSamples);
 
        fprintf('  EVM:\n');
        fprintf('    EVM peak: %0.3f%%  EVM RMS: %0.3f%%\n\n', ...
        evm.Peak,evm.RMS);
 
        fprintf('  Decoded MAC Sequence Control field contents:\n');
        fprintf('    Sequence number:%d\n',packetSeq(pktInd));
    end

% �����������
    % Update search index
    searchOffset = pktOffset+double(indNonHTData(2));
 
pktInd = pktInd+1;

% ���ظ��İ���⵽ʱ����������
    if length(unique(packetSeq))<length(packetSeq)
        break
    end  
end

packetSeq

