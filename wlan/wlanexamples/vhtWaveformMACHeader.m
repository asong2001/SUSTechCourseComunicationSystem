function mpduHeaderBits = vhtWaveformMACHeader()
% vhtWaveformMACHeader Featured example helper function
% Creates bits for a MAC header with no contents

% Copyright 2015 The MathWorks, Inc.

%#codegen

type = 2;    % Data frame type 2 (10)
subtype = 0; % Data subtype 0 (00)

% Create MPUD header
mac = struct;
frameControl = getFrameControl(type,subtype);
fields = fieldnames(frameControl);
frameControlBits = [];
for i=1:numel(fields)
    frameControlBits = [frameControlBits frameControl.(fields{i})]; %#ok<AGROW>
end
mac.FrameControl = bi2de(reshape(frameControlBits,8,[]).').'; % 2 octets
mac.Duration = [0 0];         % Duration of frame for NAV (2 octets)
mac.Address1 = [0 0 0 0 0 0]; % Destination address (6 octets) 
mac.Address2 = [0 0 0 0 0 0]; % Station address (6 octets)
mac.Address3 = [0 0 0 0 0 0]; % BSSID (6 octets)
mac.Sequence = [0 0];         % 2 octets
mac.Address4 = [0 0 0 0 0 0]; % 6 octets
mac.QoS = [0 0];              % 2 octets
mac.VHTControl = [128 0 0 0]; % 4 octets
% Convert mac header structure to bit vector
mpduHeaderBits = octetStruct2bits(mac);
end
  

% Frame control fields
function frameCtrl = getFrameControl(type,subtype)
    frameCtrl = struct;
    frameCtrl.ProtocolVersion = uint8(de2bi(0,2));
    frameCtrl.Type            = uint8(de2bi(type,2));
    frameCtrl.Subtype         = uint8(de2bi(subtype,4));
    frameCtrl.ToDS            = uint8(0);
    frameCtrl.FromDS          = uint8(0);
    frameCtrl.MoreFragments   = uint8(0);
    frameCtrl.Retry           = uint8(0);
    frameCtrl.PowerManagement = uint8(0);
    frameCtrl.MoreData        = uint8(0);
    frameCtrl.ProtectedFrame  = uint8(0);
    frameCtrl.Order           = uint8(0);
end

% Convert structure with octets to bitstream
function bits = octetStruct2bits(str)
    fnames = fieldnames(str);
    numFields = numel(fnames);
    bits = [];
    for f = 1:numFields
        octets = str.(fnames{f});
        for i=1:numel(octets) 
            bits = [bits; de2bi(double(octets(i)),8).']; %#ok<AGROW>
        end
    end
end