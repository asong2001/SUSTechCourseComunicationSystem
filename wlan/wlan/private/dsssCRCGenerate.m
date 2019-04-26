function CRC = dsssCRCGenerate(in)
%dsssCRCGenerate DSSS CRC generation
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   CRC = dsssCRCGenerate(IN)
%

%   Copyright 2015 The MathWorks, Inc.

%#codegen
    
    genPoly  =  uint16(4129); % x^16 + x^12 + x^5 + 1
    numBits  =      int8(16); % 16-bit CRC
    shiftReg = uint16(65535); % all ones
    finalXOR = uint16(65535); % all ones

    for i = 1:length(in)
        regOut = bitand(bitshift(shiftReg, 1 - numBits), uint16(1));
        regIn  = bitxor(regOut, uint16(in(i)));

        if regIn
            shiftReg = bitshift(bitxor(bitshift(genPoly, -1), shiftReg), 1);
            shiftReg = shiftReg + 1;
        else
            shiftReg = bitshift(shiftReg, 1);
        end
    end

    shiftReg = bitxor(shiftReg, finalXOR);

    CRC = coder.nullcopy(zeros(numBits, 1, 'int8'));
    for i = 1:numBits
        CRC(i) = bitand(bitshift(shiftReg, -(numBits - i)), 1);
    end

end

% [EOF]