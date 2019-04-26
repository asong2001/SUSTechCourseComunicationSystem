function refSym = helperClosestConstellationPoint(sym,cfgFormat)
%helperClosestConstellationPoint Find the closest constellation point
%
%   REFSYM = helperClosestConstellationPoint(SYM,CFGVHT) returns the
%   closest constellation point for a given symbol SYM and format
%   configuration object CFGFORMAT. Only OFDM formats are supported.
%
%   SYM is an array containing complex symbols of type double.
%
%   CFGFORMAT is the format configuration object of type <a href="matlab:help('wlanVHTConfig')">wlanVHTConfig</a>, 
%   <a href="matlab:help('wlanHTConfig')">wlanHTConfig</a>, or <a href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a>, which specifies the parameters for
%   the VHT, HT-Mixed, and Non-HT formats, respectively.
%
%   % Example: Plot the closest point for a noisy QPSK constellation
%
%     sym = awgn(1/sqrt(2)*qammod(randi([0 3],100,1),4),20); % Noisy QPSK
%     cfgVHT = wlanVHTConfig('MCS',1); % MCS1 = QPSK
%     refSym = helperClosestConstellationPoint(sym,cfgVHT);
%     figure;
%     plot(sym,'bo');
%     hold on;
%     plot(refSym,'rx');
%     legend('Symbols','Ref','Location','South');

%   Copyright 2015 The MathWorks, Inc.

validateattributes(sym,{'double'},{},mfilename,'symbol');

% Get the constellation used for the transmission
ref = helperConstellationSymbols(cfgFormat);

% Determine the closest reference constellation point for each index
refSym = zeros(size(sym));
for i = 1:numel(sym)
    [~,idx] = min(abs(sym(i)-ref));
    refSym(i) = ref(idx);
end
end