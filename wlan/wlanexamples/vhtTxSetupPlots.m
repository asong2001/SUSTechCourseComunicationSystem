function [hSF,hCon,hEVM] = vhtTxSetupPlots(cfgVHT)
%vhtTxSetupPlots Create measurement plots for
%VHTTransmitterMeasurementsExample featured example

% Copyright 2015 The MathWorks, Inc.

% Scale and poisition plots on the screen
res = get(0,'ScreenSize');
if (res(3)>1280)
    xpos = fix(res(3)*[1/2 3/4 1/4 1/4]);
    ypos = fix(res(4)*[1/16 1/2]);
    xsize = xpos(2)-xpos(1)-20;
    ysize = fix(xsize*5/6);
    repositionPlots = true;
else
    repositionPlots = false;
end

% Spectral flatness diagram
hSF = figure; 
title('Spectral Flatness, Packet: 1');
grid('on');
hSF.Visible = 'Off';
if repositionPlots
    hSF.Position = [xpos(1) ypos(2) xsize ysize];
end

% Number of spatial streams
numSS = cfgVHT.NumSpaceTimeStreams/(double(cfgVHT.STBC)+1);

% Reference constellation symbols
refConst = helperConstellationSymbols(cfgVHT);

% Constellation diagram and EVM per subcarrier plot for each spatial stream
hCon = cell(numSS,1);
hEVM = cell(numSS,1);
for i = 1:numSS
    % Constellation diagram per spatial stream
    hCon{i} = comm.ConstellationDiagram;
    hCon{i}.ReferenceConstellation = refConst;
    hCon{i}.Title = 'Equalized Data Symbols, Packet:1, Spatial Stream:1';
    
    % EVM per subcarrier per spatial stream
    hEVM{i} = figure;
    title('RMS EVM, Packet:1, Spatial Stream:1');
    grid('on');
    hEVM{i}.Visible = 'Off';
    if repositionPlots
        hCon{i}.Position = [xpos(2) ypos(2) xsize ysize];
        hEVM{i}.Position = [xpos(2) ypos(1) xsize ysize];
    end
end

end

