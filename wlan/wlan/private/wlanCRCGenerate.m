function y = wlanCRCGenerate(x)
%wlanCRCGenerate Generate CRC checksum
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
% 
%   Y = wlanCRCGenerate(X) generates CRC checksum for an input message X.
%   
%   Y is an int8 column vector of length 8 containing the checksum.
%  
%   X is an integer, binary, column vector containing the message bits.
%
%   See also wlanCRCDetect.

% Copyright 2015 The MathWorks, Inc.

%#codegen

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
