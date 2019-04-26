function vhtSigRecDisplaySIGAInfo(cfgVHTRx)
%vhtSigRecDisplaySIGAInfo Featured example helper function
%
%   Displays the recovered parameters in VHT-SIG-A

%   Copyright 2015 The MathWorks, Inc.

vhtSIGADisp = struct;
vhtSIGADisp.ChannelBandwidth = cfgVHTRx.ChannelBandwidth;
vhtSIGADisp.NumSpaceTimeStreams = cfgVHTRx.NumSpaceTimeStreams;
vhtSIGADisp.STBC = cfgVHTRx.STBC;
vhtSIGADisp.MCS = cfgVHTRx.MCS;
vhtSIGADisp.ChannelCoding = cfgVHTRx.ChannelCoding;
vhtSIGADisp.GuardInterval = cfgVHTRx.GuardInterval;
vhtSIGADisp.GroupID = cfgVHTRx.GroupID;
vhtSIGADisp.PartialAID = cfgVHTRx.PartialAID;
vhtSIGADisp.Beamforming = cfgVHTRx.Beamforming;
vhtSIGADisp.PSDULength = cfgVHTRx.PSDULength;

disp('  Decoded VHT-SIG-A contents: ');
disp(vhtSIGADisp);

end