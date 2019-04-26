function y = wlanScramble(msg,iniState)
%wlanScramble Scramble and descramble the binary input
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   Y = wlanScramble(MSG,INISTATE) scramble and descramble the binary
%   input (MSG) using a synchronous scrambler.
%
%   Y is a signed 8-bit output of size N-by-1, where N is the size of input
%   MSG.
%
%   MSG is a column vector that includes SERVICE, PSDU, tail and pad bits
%   which are scrambled with a length-127 frame-synchronous scrambler. The
%   frame-synchronous scrambler uses the generator polynomial defined in
%   IEEE(R) Std 802.11(TM)-2007. The same scrambler/descrambler is used for
%   scrambling the transmit and descramble the receive data.
%
%   The INISTATE is a row vector of size 1-by-7 of binary bits used to
%   set the initial state of the scrambler.
%   
%   Example: 
%   % Recover the input bits after scrambling and descrambling operation.
%
%   scramInit = 93;  % Scramble initialization seed
%   scramInitBits = de2bi(scramInit,7,'left-msb');
%   rng(0);
%   input = randi([0,1],1000,1);
%   scrambData = wlanScramble(input,scramInitBits);
%   descrambleData = wlanScramble(scrambData,scramInitBits);
%   isequal(input,descrambleData);
%
%   See also wlanBCCEncode.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

if isscalar(iniState) 
    iniStateVec = iniState.*ones(1,7);
else 
    iniStateVec = reshape(iniState, 1, []);
end
val = iniStateVec;
buffSize = 127;

I = zeros(1,buffSize, 'int8'); 
for d = 1:buffSize
    I(d) = xor(val(4), val(7));
    val(end:-1:2) = fliplr(val(1:end-1));
    val(1) = I(d);
end

% The scrambled sequence is handled through a circular buffer
k = 0;
y = zeros(length(msg),1, 'int8'); 
for w = 1:length(msg)  
    k = mod(k,buffSize)+1;    % Index to cycle round the scramble sequence
    y(w) = xor(msg(w),I(k));     
end

end

% [EOF]
