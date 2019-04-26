function y = wlanWindowing(x, FFTLen, wLength, guardInterval, preambleLength)                  
%wlanWindowing Window time-domain OFDM symbols
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   Y = wlanWindowing(X,FFTLEN,WLENGTH,GUARDINTERVAL,PREAMBLELENGTH)
%   returns the time-domain windowed signal for the OFDM signal. The
%   windowing function for OFDM waveform is defined in IEEE Std
%   802.11-2012.
%
%   Y is a complex Ns-by-Nt matrix containing the windowed signal, where Ns
%   is the number of time domain samples, and Nt is the number of transmit
%   antennas.
%
%   X is a complex Ns-by-Nt matrix array containing the time-domain OFDM
%   signal. 
%
%   FFTLEN is the FFT length used to modulate the OFDM symbols. The FFT
%   length depends on channel bandwidth of the input signal and should be
%   equal to 64, 128, 256 and 512 for 'CBW20', 'CBW40', 'CBW80' and
%   'CBW160' respectively. 'CBW10' and 'CBW5' also use an FFT length of 64.
%
%   WLENGTH is the windowing length in samples to apply. The WLENGTH of
%   zero samples means no windowing.
%
%   GUARDINTERVAL is the cyclic prefix for the data field within a packet,
%   specified as 'Long' or 'Short'. 
%
%   PREAMBLELENGTH is the number of samples in the preamble field of the
%   input signal.
%
%    Example:
%    %  Window a time domain waveform for an 802.11ac VHT transmission.
%
%       cfgVHT = wlanVHTConfig;       % Create format configuration
%       lstf = wlanLSTF(cfgVHT);      % Preamble generation
%       lltf = wlanLLTF(cfgVHT);
%       lsig = wlanLSIG(cfgVHT);
%       vhtsiga = wlanVHTSIGA(cfgVHT);
%       vhtstf = wlanVHTSTF(cfgVHT);
%       vhtltf = wlanVHTSTF(cfgVHT);
%       vhtsigb = wlanVHTSIGB(cfgVHT);
%       preamble = [lstf;lltf;lsig;vhtsiga;vhtstf;vhtltf;vhtsigb];
%       % Generate Data field
%       rng(0);
%       txPSDU = randi([0,1], cfgVHT.PSDULength*8,1);
%       data = wlanVHTData(txPSDU,cfgVHT);
%       txWaveform = [preamble; data];
%       % Window transmit signal
%       wLength = 16;   % Windowing samples for Transition time of 100nsec
%       sizeFFT = 256;  % FFT size for 80MHz bandwidth
%       windowedWaveform = wlanWindowing(txWaveform,sizeFFT, ...
%                      wLength,cfgVHT.GuardInterval,length(preamble));
%       % Display the spectrum
%       sa = dsp.SpectrumAnalyzer('SampleRate',80e6,'ShowLegend',true, ...
%            'Window','Flat Top','ChannelNames', ...
%            {'Transmit Waveform','Windowed waveform'});
%       % Oversample the signal by a factor of 2 before the display
%       step(sa,[resample([txWaveform;zeros(15,1)],2,1), ...
%               resample(windowedWaveform,2,1)]);
% 
%   See also wlanWaveformGenerator, wlanTGnChannel.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

validateattributes(x,{'double'},{'2d','finite','nonempty'}, ...
                           mfilename,'signal input');
if wLength == 0
    y = x;
    return;
end

% Initialize parameters
switch FFTLen
     case 128
        shortCyclicPrefix = 16;
        sampleRate = 40e6;
     case 256
        shortCyclicPrefix = 32;
        sampleRate = 80e6;
    case 512
        shortCyclicPrefix = 64;
        sampleRate = 160e6;
    otherwise % 64 for 20/10/5 MHz
        shortCyclicPrefix = 8;
        sampleRate = 20e6;  % nominal 20 MHz
end
longCyclicPrefix = FFTLen/4;
longSymLength = 4e-6*sampleRate;    % In samples, works for 5/10/20 MHz
shortSymLength = 3.6e-6*sampleRate;

% Number of transmit antennas
numTx = size(x,2);
y = complex(zeros(0,numTx));

% Standard defined windowing indices and magnitude for Long symbols
[indLong, magLong] = windowingEquation(wLength,longSymLength);
ofdmSymLen = FFTLen+longCyclicPrefix;

if strcmp(guardInterval,'Long')
  
    %Check the wLength to be less than or equal to twice the CP
    coder.internal.errorIf (strcmpi(guardInterval,'Long') && ...
                       (wLength > 2*longCyclicPrefix), ...
                       'wlan:wlanWindowing:InvalidWindowLength');
                   
    %Windowing L-STF and L-LTF fields
    nonHTtrainingSym = x(1:4*(ofdmSymLen),:); 
    remainingSym = x(4*(ofdmSymLen)+1:end,:);
    
    nonHTSymExtended = windowLSTFandLLTF(nonHTtrainingSym,FFTLen, ...
                         longCyclicPrefix,longSymLength,wLength,numTx);
    
    % Process remaining symbols
    numSym =  length(remainingSym)/(ofdmSymLen);
     
    % Check input length
    coder.internal.errorIf (mod(numSym,1)~=0, ...
                             'wlan:wlanWindowing:IncorrectInputLength');
   
    % Reshape by SymLength-NumSym-NumTx
    remainingSym = reshape(remainingSym,(ofdmSymLen), ...
                                numSym,numTx);
     
    % Extend symbol length before windowing
    remainingSymExt = [remainingSym(end-(abs(indLong(1))+FFTLen/2)+1:...
           end-FFTLen/2,:,:);...
           remainingSym;remainingSym(FFTLen/2+(1:wLength/2),:,:)];
       
    % Apply windowing on the extended symbol portion 
    remaingSymExtWin =  bsxfun(@times,remainingSymExt,magLong);
     
    nonHTSymOut = windowSym(nonHTSymExtended,wLength);
    suffixNonHTSym = nonHTSymOut(end-(wLength)+2:end,:);
    remaingSymOut = windowSym(remaingSymExtWin,wLength);
    prefixSym = remaingSymOut(1:wLength-1,:);
    y = [nonHTSymOut(1:end-wLength+1,:); ...
        (suffixNonHTSym+prefixSym);  ....
         remaingSymOut(wLength:end,:)];
    
elseif strcmp(guardInterval, 'Short')
    
    % Check the wLength to be less than or equal to twice the CP              
    coder.internal.errorIf (strcmpi(guardInterval,'Short') && ...
                       (wLength > 2*shortCyclicPrefix), ...
                       'wlan:wlanWindowing:InvalidWindowLength')
    % Extract data field
    data = x(preambleLength+1: end,:);
    
    %Windowing L-STF and L-LTF fields
    firstTwoSym = x(1:4*ofdmSymLen,:);
    nonHTSymExtended = windowLSTFandLLTF(firstTwoSym,FFTLen, ...
                         longCyclicPrefix,longSymLength,wLength,numTx);
    
    % Reshape preamble (less L-LTF and L-STF) by SymLength-NumSym-NumTx
    preambleLessFirstTwoFields = x(4*ofdmSymLen+1:preambleLength,:);
    numSym = size(preambleLessFirstTwoFields,1)/longSymLength; 
    
    preambleSym = reshape(preambleLessFirstTwoFields,...
                  FFTLen+longCyclicPrefix, ...
                  numSym,numTx);
    
    % Extend symbol length
    preambleSymExt = [preambleSym(end-(abs(indLong(1))+FFTLen/2) ...
                    +1:end-FFTLen/2,:,:);preambleSym; ...
                    preambleSym(FFTLen/2+(1:wLength/2),:,:)];
    
    % Apply windowing to all preamble symbols
    remaingSymExtWin = bsxfun(@times,preambleSymExt,magLong);

    nonHTSymOut = windowSym(nonHTSymExtended,wLength);
    suffixNonHTSym = nonHTSymOut(end-(wLength)+2:end,:);
     
    preambleSym = windowSym(remaingSymExtWin,wLength);
    prefixPreambleSym = preambleSym(1:wLength-1,:);
    suffixPreambleSym = preambleSym(end-(wLength)+2:end,:);

    % Number of data symbols
    numDataSym = length(data)/(FFTLen+shortCyclicPrefix);
    
    % Check input length
    coder.internal.errorIf (mod(numDataSym,1)~= 0, ...
                            'wlan:wlanWindowing:IncorrectInputLength');
    
    if isempty(data) % NDP
        % Only the preamble fields
        y = [nonHTSymOut(1:end-wLength+1,:); ...
        (suffixNonHTSym+prefixPreambleSym); ...
        preambleSym(wLength:end-wLength+1,:); ...
        suffixPreambleSym];               
    else
        % Window the data symbols
        % Reshape data by SymLength x  NumSym x NumTx
        dataSym = reshape(data,FFTLen+shortCyclicPrefix,numDataSym,numTx);
    
        % Standard defined windowing equations
        [index, magnitude] = windowingEquation(wLength,shortSymLength);

        % Extend symbol length
        dataSymExt = [dataSym(end-(abs(index(1))+FFTLen/2)+1: ...
                 end-FFTLen/2,:,:);dataSym;dataSym(FFTLen/2+(1:wLength/2),:,:)];

        % Apply windowing to all symbols
        dataSymExtWin = bsxfun(@times,dataSymExt,magnitude);
        dataSym = windowSym(dataSymExtWin,wLength);
        prefixdataSym = dataSym(1:wLength-1,:);

        y = [nonHTSymOut(1:end-wLength+1,:); ...
            (suffixNonHTSym+prefixPreambleSym); ...
            preambleSym(wLength:end-wLength+1,:); ...
            (suffixPreambleSym+prefixdataSym);...
            dataSym(wLength:end,:)];               
    end
    
end

end

function out = windowSym(inputSym, wlength)
% This function is performing the overlap and add operation on the extended
% window potion.

% Set the windowing length equal to the preIndx window length.
W = wlength-(wlength>0);

firstRegion = inputSym(1:end-W,1,:);
overlapRegion = inputSym(end-W+1:end,1:end-1,:)+inputSym(1:W,2:end,:);
middleRegion = [overlapRegion; inputSym(W+1:end-W,2:end,:)];
lastRegion = inputSym(end-W+1:end,end,:);
out = [firstRegion;reshape(middleRegion,[],1,size(inputSym,3));lastRegion];

end

function [windowIdx,windowMag] = windowingEquation(TTR, symLength)
% The windowing function for the OFDM symbols is defined in IEEE Std
% 802.11-2012, Section 18.3.2.5, Equation 18-4.
    preIdx = -TTR/2+1:TTR/2-1;
    midIdx = TTR/2:(symLength-TTR/2)-1;
    postIdx = symLength-TTR/2:(symLength+TTR/2)-1;
    windowIdx = [preIdx midIdx postIdx].'; 
    preMag = sin(pi/2*(0.5+(preIdx)/(TTR))).^2;
    midMag = ones(1,length(midIdx)); 
    postMag = sin(pi/2*(0.5-((postIdx)-symLength)/TTR)).^2; 
    windowMag = [preMag midMag postMag].';

end

function out = windowLSTFandLLTF(preamble,sizeFFT,longCyclicPrefix, ...
                            longSymLength,wLength,numTx)
    
    ofdmSymLen = sizeFFT+longCyclicPrefix;
    % Process first two preamble symbols
    nonHTtrainingSym = preamble(1:4*(ofdmSymLen),:);
    
    % Standard defined windowing equations
    [index,magnitude] = windowingEquation(wLength,2*longSymLength);

    % Number of initial symbols
    numSym =  length(nonHTtrainingSym)/(2*(ofdmSymLen));
   
    % Reshape by SymLength movingAvgLen  NumSym movingAvgLen NumTx
    nonHTtrainingSym = reshape(nonHTtrainingSym,2* ...
                                    (ofdmSymLen),numSym,numTx);
    
    % Extend symbol length
    initialSymExtended = [nonHTtrainingSym(end-(abs(index(1))+sizeFFT/2) ...
                    + 1:end-sizeFFT/2,:,:);nonHTtrainingSym; ...
                    nonHTtrainingSym(sizeFFT/2+(1:wLength/2),:,:)];
    
    % Apply windowing to first two preamble symbols
    out = bsxfun(@times,initialSymExtended,magnitude);
end

% [EOF]
