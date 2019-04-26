function [psdu,cfgVHT] = vhtWaveformGeneratePSDU(mpdu,cfgVHT)
% vhtWaveformGeneratePSDU Featured example helper function
% Pads and delimits an MPDU to form a PSDU according to IEEE Std
% 802.11ac-2013 Section 9.12.6

% Copyright 2015 The MathWorks, Inc.

%#codegen

% Only support VHT single user configuration
validateattributes(cfgVHT,{'wlanVHTConfig'},{'scalar'},mfilename, ...
    'format configuration object');
if cfgVHT.NumUsers>1
    error('Only single user operation supported');
end
validateattributes(mpdu,{'numeric'},{'nonempty','column'},mfilename,'MPDU');

% Set A-MPDU length pre EOF padding
numMPDUDelimiterOctets = 4;
cfgVHT.APEPLength = numMPDUDelimiterOctets+numel(mpdu)/8; 

% Calculate required padding to PSDULength as per IEEE Std
% 802.11ac-2013 Section 9.12.6
ampduLength = cfgVHT.APEPLength;
numPadOctetsPreEOF = 0;
while (ampduLength<cfgVHT.PSDULength(1))&&(mod(ampduLength,4)~=0)
    numPadOctetsPreEOF = numPadOctetsPreEOF+1;
    ampduLength = ampduLength+1;
end
numEOFDelimiters = 0;
while (ampduLength+4)<=cfgVHT.PSDULength(1)
    numEOFDelimiters = numEOFDelimiters+1;
    ampduLength = ampduLength+4;
end
numPadOctetsPostEOF = 0;
while ampduLength<cfgVHT.PSDULength(1)
    numPadOctetsPostEOF = numPadOctetsPostEOF+1;
    ampduLength = ampduLength+1;
end

% Generate delimiter for A-MPDU. Depending on the EOF padding it may be
% the only subfarme and therefore a VHT single MPDU. In this case set
% EOF field in delimiter to 1.
if numEOFDelimiters == 0
    eof = 1; % VHT single MPDU
else
    eof = 0;
end
mpduDelimiter = delimiter(mpdu,eof); % Generate MPDU delimiter
eofDelimiter = delimiter([],1);      % Generate EOF delimiter
eofDelimiter = repmat(eofDelimiter,numEOFDelimiters,1);
prePad = zeros(numPadOctetsPreEOF*8,1);
postPad = zeros(numPadOctetsPostEOF*8,1);

% Form PSDU with A-MPDU, EOF-delimiters and padding
psdu = [mpduDelimiter; mpdu; prePad; eofDelimiter; postPad];

end

% Creates an MPDU or EOF delimiter 
function bits = delimiter(mpdu,eof)
    % IEEE Std 802.11ac-2013 Section 8.6.1    
    mpduLength = de2bi(numel(mpdu)/8,14).'; % 14 bits in length field
    preCRCBits = [eof; 0; mpduLength]; % 0 is reserved bit
        
    % CRC generation
    genMPDU = comm.CRCGenerator([8 2 1 0], ...
        'InitialConditions',1,'DirectMethod',true,'FinalXOR',1);
    postCRCBits = step(genMPDU,preCRCBits); % Append CRC
    
    % Append signature
    signature = de2bi(hex2dec('4E'),8).';
    bits = [postCRCBits; signature];
end
