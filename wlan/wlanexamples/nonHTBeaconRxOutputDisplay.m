function nonHTBeaconRxOutputDisplay(mpdu)
% Helper function for Non-HT beacon receiver example. 

% Copyright 2015 The MathWorks, Inc.

%#codegen

body = mpdu.FrameBody;

% Determine if AP is in HT mode (11n/11ac)
HTCapable = false;
for idx = 1:length(body.InfoElements)
   if body.InfoElements(idx).ID == 45
      HTCapable = true; 
      break;
   end 
end

% Get SSID
SSID = char(body.InfoElements(1).Value(1:body.InfoElements(1).Length));

% Print out information
if strcmpi(SSID,'') || isempty(SSID)
    disp('SSID: (Blank)');
else
    disp(['SSID: ', SSID]);
end

if HTCapable
    disp('HT Capable');
end

disp('------------------------');

end

% [EOF]