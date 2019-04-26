function vhtTxEVMConstellationPlots(eqSym,subcarrierEVM,cfg,pktNum,hCon,hEVM)  
% VHTTransmitterMeasurements Featured example helper function
% Plots EVM per subcarrier and constellation.

% Copyright 2015 The MathWorks, Inc.

% Function for VHT only
validateattributes(cfg,{'wlanVHTConfig'},{'scalar'},mfilename, ...
    'format configuration object');

numSpatialStreams = cfg.NumSpaceTimeStreams/(double(cfg.STBC)+1);

Nfft = helperFFTLength(cfg);
dataIndices = helperSubcarrierIndices(cfg,'VHT');
for i = 1:numSpatialStreams
    % Plot equalized data symbols
    tmp = eqSym(:,:,i);
    step(hCon{i},tmp(:))
    hCon{i}.Title = ['Equalized Data Symbols, Packet:' num2str(pktNum) ...
        ', Spatial Stream:' num2str(i)];

    figure(hEVM{i});
    hEVMAx = gca;
    
    % Get the previous upper YLimit and if the maxium EVM in the current
    % packet is greater then increase the limit
    prevLim = hEVMAx.YLim(2);
    currentLim = ceil(max(subcarrierEVM(:,:,i)));
    if currentLim>prevLim
        lim = [0 currentLim];
    else
        lim = [0 prevLim];
    end
    
    % Plot the RMS EVM accross subcarriers
    plotInd = dataIndices-Nfft/2-1;
    plot(hEVMAx,plotInd,subcarrierEVM(:,:,i),'bo');
    ylim(hEVMAx,lim);
    title(hEVMAx,['RMS EVM, Packet:' num2str(pktNum) ...
        ', Spatial Stream:' num2str(i)]);
    ylabel(hEVMAx,'EVM (%)');
    xlabel(hEVMAx,'Subcarrier Index');
    grid(hEVMAx,'on');
end
drawnow;
end