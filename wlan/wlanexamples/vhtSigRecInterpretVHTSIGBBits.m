function [length,crcBits] = vhtSigRecInterpretVHTSIGBBits(sigBBits,cfgVHT)
% vhtSigRecInterpretVHTSIGBBits interprets VHT-SIG-B bits
%
%   [LENGTH,CRCBITS] = vhtSigRecInterpretVHTSIGBBits(BITS,CFGVHT) returns
%   the decoded APEP length and the reference CRC bits given recovered
%   VHT-SIG-B bits, BITS, and a VHT configuration object.

%   Copyright 2015 The MathWorks, Inc.

switch cfgVHT.ChannelBandwidth
    case 'CBW20' % 26
        length = bi2de(sigBBits(1:17).')*4;
        crcBits = sigBBits(1:20); % VHT-SIG-B excluding crc
    case 'CBW40' % 27
        length = bi2de(sigBBits(1:19).')*4;
        crcBits = sigBBits(1:21); % VHT-SIG-B excluding crc
    otherwise    % 29 for {'CBW80', 'CBW80+80', 'CBW160'}
        length = bi2de(sigBBits(1:21).')*4;
        crcBits = sigBBits(1:23); % VHT-SIG-B excluding crc
end
crcBits = generateCRC(crcBits);

end

% Generate a CRC for input bits x.
function y = generateCRC(x)

% Section 20.3.9.4.4 in IEEE Std 802.11-2012
genPoly  = uint8(7);   
numBits  = int8(8);
shiftReg = uint8(255); 
finalXOR = uint8(255);

for i = 1:length(x)
    regOut = bitand(bitshift(shiftReg, 1 - numBits), uint8(1));
    regIn  = bitxor(regOut, uint8(x(i)));

    if regIn
        shiftReg = bitshift(bitxor(bitshift(genPoly, -1), shiftReg), 1);
        shiftReg = shiftReg + 1;
    else
        shiftReg = bitshift(shiftReg, 1);
    end
end  

shiftReg = bitxor(shiftReg, finalXOR);         

y = coder.nullcopy(zeros(numBits, 1, 'int8'));
for i = 1:numBits
    y(i) = bitand(bitshift(shiftReg, -(numBits - i)), 1);
end

end