function [SpectrumAnalyzer,TimeScope,ConstellationDiagram] = ...
    vhtSigRecSetupPlots(sr)
%vhtSigRecSetupPlots Featured example helper function
%
%   Creates and configures analysis objects

%   Copyright 2015 The MathWorks, Inc.

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

% Create spectrum analyzer to plot spectrum of packet
SpectrumAnalyzer = dsp.SpectrumAnalyzer;
SpectrumAnalyzer.Name = 'Detected packet signal spectrum';
SpectrumAnalyzer.SampleRate = sr;
SpectrumAnalyzer.ShowLegend = true;
SpectrumAnalyzer.RBWSource = 'Property';
SpectrumAnalyzer.RBW = 100e3;
SpectrumAnalyzer.SpectralAverages = 3;
if repositionPlots
    SpectrumAnalyzer.Position = [xpos(1) ypos(2) xsize ysize];
end

% Create time scope to plot detected packet
TimeScope = dsp.TimeScope;
TimeScope.Name = 'Detected packet';
TimeScope.YLabel = 'Magnitude';
TimeScope.PlotType = 'Line';
TimeScope.ReduceUpdates = false;
TimeScope.ShowGrid = true;
TimeScope.ShowLegend = true;
TimeScope.SampleRate = sr;
if repositionPlots
    TimeScope.Position = [xpos(2) ypos(2) xsize ysize];
end

% Create a constellation diagram for each potential spatial stream
NssMax = 8; % Maximum number of spatial streams
ConstellationDiagram = cell(NssMax,1);
for iss = 1:NssMax
    ConstellationDiagram{iss} = comm.ConstellationDiagram;
    ConstellationDiagram{iss}.ShowReferenceConstellation = true;
    ConstellationDiagram{iss}.ShowGrid = true;
    ConstellationDiagram{iss}.Name = 'Equalized data symbols';
    ConstellationDiagram{iss}.Title = ['Spatial stream ' num2str(iss)];
    if (repositionPlots)
        if mod(iss,2)
            ConstellationDiagram{iss}.Position = [xpos(1) ypos(1) xsize ysize];
        else
            ConstellationDiagram{iss}.Position = [xpos(2) ypos(1) xsize ysize];
        end
    end
end

end