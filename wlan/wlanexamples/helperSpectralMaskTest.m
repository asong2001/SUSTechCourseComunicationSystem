function helperSpectralMaskTest(x, fs, osr, varargin)
%helperSpectralMaskTest Featured example helper function
%   Plots the power spectral density (PSD) and overlays WLAN PSD limits to
%   check if spectral emissions are within specified levels.

%   Copyright 2015 The MathWorks, Inc.

narginchk(3,5);

cbwMHz = fs/1e6;    % Channel bandwidth in MHz
if nargin>3
    % Must be same size vectors
    dBrLimits = varargin{1};
    fLimits = varargin{2};
else % Default
    %  From IEEE Std 802.11ac-2013 Section 22.3.18.1
    dBrLimits = [-40 -40 -28 -20 0 0 -20 -28 -40 -40];
    fLimits = [-Inf -1.5*cbwMHz -cbwMHz -(cbwMHz/2+1) -(cbwMHz/2-1) ...
        (cbwMHz/2-1) (cbwMHz/2+1) cbwMHz 1.5*cbwMHz Inf];
end

% Handle edges when spectral mask is flat
if (-cbwMHz*osr/2)<fLimits(2) && (cbwMHz*osr/2)>fLimits(end-1)
    fLimits(1) = -fs/1e6*osr/2;
    fLimits(end) = fs/1e6*osr/2;
end

% Only plot or compare at frequencies applicable to the sampling rate
plotIdx = (fLimits>=((-fs*osr/2)/(1e6)))&(fLimits<=((fs*osr/2)/(1e6)));

rbw = 100e3;        % Resolution bandwidth
vbw = 30e3;         % Video bandwidth
N = floor(rbw/vbw); % number of spectral averages

% Calculate segment length based on rbw, sampling rate and window.
segLen = calcSegmentLength(rbw, fs*osr, 'Hann');
if length(x) < segLen
    error('Not enough input data to measure spectral mask.');
end
numSegments = floor(length(x(:,1))/segLen);

% Switch to dsp.SpectrumAnalyzer once mask levels are enabled.
SE = dsp.SpectrumEstimator('SampleRate', fs*osr,...
    'SpectrumType','Power density',...
    'PowerUnits', 'dBm',...
    'FrequencyRange','centered',...
    'SpectralAverages', N ,...
    'FFTLengthSource', 'Property', 'FFTLength', segLen);

psdAv = zeros(segLen, numSegments);
for txIdx = 1:size(x,2)
    figure;
    for idx = 1:numSegments
        psdAv(:,idx) = step(SE, x((idx-1)*segLen+(1:segLen),txIdx));
    end    
    freqs = getFrequencyVector(SE);

    % Get the dBr limit for each frequency bin
    dBrBinLimit = interp1(fLimits(plotIdx), dBrLimits(plotIdx), ...
                          freqs/1e6);

    % Upper Y limit for plotting
    upLim = -(floor(abs(max(psdAv(:)))/10))*10;
    % Lower Y limit for plotting (limit to at most 90dB from upper limit)
    lowLim = max(-(ceil(abs(min(psdAv(:)))/10))*10,upLim-90);

    for idx = 1:numSegments
        % Plot each PSD computation
        plot(freqs/1e6, psdAv(:,idx), 'b');
        hold on
        plot(fLimits(plotIdx), dBrLimits(plotIdx)+max(psdAv(:,idx)), 'o-r');
        ylim([lowLim upLim]);
        legend({['Transmit antenna ' num2str(txIdx)], 'Relative mask'}, ...
            'Location','SouthEast');
        grid on;
        title(['Spectral Mask, Antenna ' num2str(txIdx)]);
        xlabel('Frequency (MHz)');
        ylabel('PSD (dBm)');
        hold off;
        drawnow;
        
        % Calculate if the mask is exceeded in any bin
        failBins = find((psdAv(:,idx)-max(psdAv(:,idx)))>dBrBinLimit);
        if any(failBins)
            disp('   Spectrum mask failed.');
            return; % Do not update any more
        end
    end    
end
disp('   Spectrum mask passed.');

end

%--------------------------------------------------------------------------
function [pSegmentLength,ActualRBW] = calcSegmentLength(desiredRBW,actualSampleRate,win)
% Get segment length corresponding to specified RBW.

% Segment length will depend on ENBW (which in turns depends on segment
% length). Thus, an initial ENBW is obtained using a segment length of
% 1000
ENBW = getENBW(1000, win);
% Compute segment length
segLen = ceil(ENBW*actualSampleRate/desiredRBW);
% Iterate over segment length to minimize error between requested RBW and
% actual RBW
count = 1;
segLenVect = segLen;
while (count<100)  % protect against very long convergence
    new_segLen = ceil( getENBW(ceil(segLen), win)*actualSampleRate/desiredRBW);
    err = abs(new_segLen - segLen);
    if (err==0)  % we have converged
        segLen = new_segLen;
        break
    end
    if ~any(segLenVect==new_segLen)
        segLenVect = [segLenVect, new_segLen]; %#ok<AGROW>
        segLen = new_segLen;
        count = count + 1;
    else
        % We hit a previously computed segment length. The sequence
        % will repeat. Break out and select the segment length that
        % minimizes the error
        L = length(segLenVect);
        computed_RBW = zeros(L, 1);
        for ind = 1:L
            % Get RBW corresponding to segLenVect(ind)
            computed_RBW(ind) = getENBW(segLenVect(ind), win)*actualSampleRate/segLenVect(ind);
        end
        % Select segment length that minimizes absolute error between
        % actual and desired RBW:
        RBWErr = abs(desiredRBW - computed_RBW);
        [~,ind_min] = min(RBWErr);
        segLen = segLenVect(ind_min);
        break
    end
end

pSegmentLength = segLen;
ActualRBW = getENBW(segLen, win)*actualSampleRate/segLen;

end

%--------------------------------------------------------------------------
function [ENBW, winFcn] = getENBW(L, Win)
% Get window parameters based on a segment length L
if isempty(L)
    % Segment length not computed yet
    L = 1000;
end
switch Win
    case 'Hann'
        w = hann(L, 'periodic');
        ENBW = (sum( w .^ 2 )/sum( w )^2)*L;
        winFcn = @(x) hann(x, 'periodic');
end
end

% [EOF]