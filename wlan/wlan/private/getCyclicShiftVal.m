function out = getCyclicShiftVal(formatType, Ntx, chBW)
%getCyclicShiftVal Get cyclic shift values
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   cSh = getCyclicShiftVal(formatType, Ntx, chBW)
%       where formatType = 'OFDM' or 'VHT'
%                   chBW = 20/40/80/160 in MHz
%
%   Output shift values are in number of samples (and not time).
%
%   As per IEEE Std 802.11ac-2013, Table 22-10, for L-STF, L-LTF, 
%       L-SIG, VHT-SIG-A, HT-SIG fields
%   As per IEEE Std 802.11ac-2013, Table 22-11, for VHT-STF, VHT-LTF,
%       VHT-SIG-B, HT-STF, HT-LTF and DATA fields.
%
%   See also wlanCyclicShift.

%   Copyright 2015 The MathWorks, Inc.

%#codegen
if strcmp(formatType, 'VHT')
    switch Ntx
        case 2
            cShift = [0; -400];
        case 3
            cShift = [0; -400; -200];
        case 4
            cShift = [0; -400; -200; -600];
        case 5
            cShift = [0; -400; -200; -600; -350];
        case 6
            cShift = [0; -400; -200; -600; -350; -650];
        case 7
            cShift = [0; -400; -200; -600; -350; -650; -100];
        case 8
            cShift = [0; -400; -200; -600; -350; -650; -100; -750];
        otherwise
            cShift = 0;
    end
else 
    switch Ntx
        case 2
            cShift = [0; -200];
        case 3
            cShift = [0; -100; -200];
        case 4
            cShift = [0; -50; -100; -150];
        case 5
            cShift = [0; -175; -25; -50; -75];
        case 6
            cShift = [0; -200; -25; -150; -175; -125];
        case 7
            cShift = [0; -200; -150; -25; -175; -75; -50];
        case 8
            cShift = [0; -175; -150; -125; -25; -100; -50; -200];
        otherwise
            cShift = 0;
    end
end

% Cyclic shift delay in number of samples. 
out = cShift*1e-9/(1/(chBW*1e6));

end

% [EOF]