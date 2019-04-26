function delay = helperCyclicShiftDelay(Nsts,fs,fieldType)
% helperCyclicShiftDelay cyclic shifts to apply to streams in samples
%
%   DELAY = helperCyclicShiftDelay(NUMSTS,FS,FIELDTYPE) returns the delay
%   in samples to apply to each stream in samples as a N-by-1 column
%   vector, where N is the number of space-time streams or transmit
%   antennas.
%
%   NUMSTS is the number of space-time streams or transmit antennas if
%   the cyclic shift is for a HT/VHT or a legacy field.
%
%   FS is the baseband sampling rate in Hertz.
%
%   FIELDTYPE is a string and must be one of 'Legacy','HT','VHT'.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

validatestring(fieldType,{'VHT','HT','Legacy'}, ...
    mfilename,'field type');
validateattributes(Nsts,{'numeric'},{'scalar','>',0,'<',9}, ...
    mfilename,'number of streams');
validateattributes(fs,{'numeric'},{'scalar'}, ...
    mfilename,'baseband sampling rate');

if any(strcmp(fieldType,{'VHT','HT'}))
    % VHT/HT portion; IEEE Std 802.11ac-2013 Table 22-11 (time in ns)
    switch Nsts
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
        otherwise % 1
            cShift = 0;
    end
else
    % NonHT portion; IEEE Std 802.11ac-2013 Table 22-10 (time in ns)
    switch Nsts
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
        otherwise % 1
            cShift = 0;
    end
end

% Cyclic shift delay in number of samples 
delay = cShift*1e-9*fs;

end