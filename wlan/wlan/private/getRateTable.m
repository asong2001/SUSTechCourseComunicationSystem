function mcsTable = getRateTable(cfgFormat)
%getRateTable Select the Rate parameters for format
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%%
%   mcsTable = getRateTable(CFGFORMAT) returns the modulation and coding
%   parameters for the format configuration object CFGFORMAT. 

%   Copyright 2015 The MathWorks, Inc.

% References:
% [1] IEEE Std 802.11ac - 2013.
% [2] "Next Generation Wireless LANs", E. Perahia & R. Stacey, Cambridge
% University Press, 2013. Section 7.4.
% [3] IEEE Std 802.11 - 2012.

%#codegen

if isa(cfgFormat, 'wlanVHTConfig') % VHT format
    
    switch cfgFormat.ChannelBandwidth
        case 'CBW20'
            NSD = 52;
        case 'CBW40'
            NSD = 108;
        case {'CBW80', 'CBW80+80'}
            NSD = 234;
        otherwise % case 'CBW160'
            NSD = 468;
    end
    
    numUsers = cfgFormat.NumUsers;
    numSS = cfgFormat.NumSpaceTimeStreams / (((numUsers == 1) && cfgFormat.STBC) + 1);
    MCS = repmat(cfgFormat.MCS, 1, numUsers/length(cfgFormat.MCS));
    
    [rate, Nbpscs, Ncbps, Ndbps, Nes] = deal(zeros(1, numUsers));
    for u = 1:numUsers
        [rate(u), Nbpscs(u), Ncbps(u), Ndbps(u), Nes(u)] = ...
            getMCSTableForOneUser(MCS(u), NSD, numSS(u));
    end

    mcsTable = struct( ...
        'Rate',       rate, ...
        'NBPSCS',     Nbpscs, ...
        'NSD',        NSD, ...
        'NCBPS',      Ncbps, ...
        'NDBPS',      Ndbps, ...
        'NES',        Nes, ...
        'Nss',        numSS);
    
elseif isa(cfgFormat, 'wlanHTConfig') % HT-Mixed format

    switch cfgFormat.ChannelBandwidth
        case 'CBW20'
            NSD = 52;
        otherwise % CBW40
            NSD = 108;
    end
    
    mcsTable = getHTMCSTable(cfgFormat.MCS, NSD);
    
elseif isa(cfgFormat, 'wlanNonHTConfig') % Non-HT format, for OFDM
    
    mcsTable = getNonHTMCSTable(cfgFormat.MCS);
    
else % Add other formats

    mcsTable = [];
end

%-------------------------------------------------------------------------
function mcsTable = getNonHTMCSTable(mcs)
% Supports data rate values in set of {6, 9, 12, 18, 24, 36, 48, 54} Mbps
% as a fcn of the MCS values.

Nsd = 48;       % Data subcarriers
switch mcs
  case 0 % 6 Mbps
    Nbpscs = 1;  % 'BPSK'
    rate = 1/2;
  case 1 % 9 Mbps
    Nbpscs = 1; 
    rate   = 3/4;
  case 2 % 12 Mbps
    Nbpscs = 2;  % QPSK
    rate   = 1/2;
  case 3 % 18 Mbps
    Nbpscs = 2; 
    rate   = 3/4;
  case 4 % 24 Mbps
    Nbpscs = 4;  % 16QAM 
    rate   = 1/2;
  case 5 % 36 Mbps
    Nbpscs = 4;  
    rate   = 3/4;
  case 6  % 48 Mbps
    Nbpscs = 6;  % '64QAM'
    rate   = 2/3;
  otherwise % 7 => 54 Mbps
    Nbpscs = 6;
    rate   = 3/4;
end    

Ncbps = Nsd * Nbpscs;
Ndbps = Ncbps * rate;  

mcsTable = struct( ...
    'Rate',       rate, ...
    'NBPSCS',     Nbpscs, ...
    'NSD',        Nsd, ...
    'NCBPS',      Ncbps, ...
    'NDBPS',      Ndbps);

%-------------------------------------------------------------------------
function mcsTable = getHTMCSTable(MCS, Nsd)
% Supports MCS values only in the range 0-31 for now.

% Strip Nss from MCS to get MC and Nss
Nss = floor(MCS/8)+1;
mc = rem(MCS, 8);

switch mc
  case 0
    Nbpscs = 1; % 'BPSK'
    rate   = 1/2;
  case 1
    Nbpscs = 2; % 'QPSK'
    rate   = 1/2;
  case 2
    Nbpscs = 2; 
    rate   = 3/4;
  case 3
    Nbpscs = 4; % '16QAM'
    rate   = 1/2;
  case 4
    Nbpscs = 4; 
    rate   = 3/4;
  case 5
    Nbpscs = 6; % '64QAM'
    rate   = 2/3;
  case 6
    Nbpscs = 6; 
    rate   = 3/4;
  otherwise % MCS == 7
    Nbpscs = 6;
    rate   = 5/6;
end    

Ncbps = Nsd * Nbpscs * Nss;
Ndbps = Ncbps * rate;  

% Any Ndbps>1200 => Nes=2, for 300 Mbps per encoder
%   Confirmed with Tables 20-30 to 20-44.
Nes = ceil(Ndbps/(4*300));  

mcsTable = struct( ...
    'Rate',       rate, ...
    'NBPSCS',     Nbpscs, ...
    'NSD',        Nsd, ...
    'NCBPS',      Ncbps, ...
    'NDBPS',      Ndbps, ...
    'NES',        Nes, ...
    'Nss',        Nss);

%-------------------------------------------------------------------------
function [rate, Nbpscs, Ncbps, Ndbps, Nes] = getMCSTableForOneUser(MCS, Nsd, Nss)

switch MCS
  case 0
    Nbpscs = 1; % 'BPSK'
    rate   = 1/2;
  case 1
    Nbpscs = 2; % 'QPSK'
    rate   = 1/2;
  case 2
    Nbpscs = 2; 
    rate   = 3/4;
  case 3
    Nbpscs = 4; % '16QAM'
    rate   = 1/2;
  case 4
    Nbpscs = 4; 
    rate   = 3/4;
  case 5
    Nbpscs = 6; % '64QAM'
    rate   = 2/3;
  case 6
    Nbpscs = 6; 
    rate   = 3/4;
  case 7
    Nbpscs = 6;
    rate   = 5/6;
  case 8
    Nbpscs = 8; % 256QAM
    rate   = 3/4;
  otherwise % MCS == 9
    Nbpscs = 8;
    rate   = 5/6;
end    

Ncbps = Nsd * Nbpscs * Nss;
Ndbps = Ncbps * rate;  

% Handle exceptions to Nes generic rule - Table 7.13 [2].
%   For each case listed, work off the Ndbps value and create a look-up
%   table for the Nes value.
%   Only 9360 has a valid value from the generic rule also, 
%   all others are exceptions
NdbpsVec = [2457 8190 9828 9360 14040 9828 16380 19656 21840 14976 22464];
expNes =   [   3    6    6    6     8    6     9    12    12     8    12];

exceptIdx = find(Ndbps == NdbpsVec);
if ~isempty(exceptIdx)
    if (Ndbps == 9360) && (Nss == 5) % One valid case for 160, 80+80
        Nes = 5;
    else  % Two exception cases
        Nes = expNes(exceptIdx(1));
    end
else  % Generic rule: 3.6*600 - for a net 600Mbps per encoder
    Nes = ceil(Ndbps/2160);
end

% [EOF]