function phi = dsssCCKModulate(scrambledPSDU,dataRate,refPhase)
%dsssCCKModulate DSSS CCK modulation
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   PHI = dsssCCKModulate(SCRAMBLEDPSDU,DATARATE,REFPHASE)
%

%   Copyright 2015 The MathWorks, Inc.

%#codegen

    % Only valid for 5.5Mbps and 11Mbps data rates
    coder.internal.errorIf(~any(strcmpi(dataRate,{'5.5Mbps','11Mbps'})), ...
        'wlan:dsssCCKModulate:InvalidDataRate');

    % Establish cckOrder, the number of data bits per CCK symbol
    if (strcmpi(dataRate,'5.5Mbps'))
        % Clause 17.4.6.6.3 CCK 5.5 Mb/s modulation
        cckOrder = 4;
    else
        % Clause 17.4.6.6.4 CCK 11 Mb/s modulation
        cckOrder = 8;
    end

    % Establish cckSymbols, the number of CCK symbols in the PSDU
    cckSymbols = (length(scrambledPSDU)/cckOrder);
    
    % Reshape the data bits to put modulation symbols in rows
    d = reshape(scrambledPSDU,cckOrder,cckSymbols).';

    % Calculate phi1
    phi1 = zeros(cckSymbols,1);
    % Initialize phi1 with reference phase from final header symbol
    phi1(1) = refPhase;
    % Extra 180 degree (pi) rotation on odd-numbered symbols
    phi1(2:2:end) = pi;
    % DQPSK modulation of dibit (d0,d1)
    dqpsk_in = [bit(d,0).'; bit(d,1).'];
    phi1 = mod(cumsum(phi1) + angle(dsssPSKModulate(dqpsk_in(:),'2Mbps')), 2*pi);
    
    % Calculate phi2, phi3 and phi4
    if (strcmpi(dataRate,'5.5Mbps'))
        % Clause 17.4.6.6.3 CCK 5.5 Mb/s modulation
        phi2 = (double(bit(d,2)) * pi) + pi/2;
        phi3 = zeros(cckSymbols,1);
        phi4 = double(bit(d,3)) * pi;
    else
        % Clause 17.4.6.6.4 CCK 11 Mb/s modulation
        phi2 = qpsk(bit(d,2),bit(d,3));
        phi3 = qpsk(bit(d,4),bit(d,5));
        phi4 = qpsk(bit(d,6),bit(d,7));
    end

    % Create overall set of phases the CCK symbols
    phi = [phi1 phi2 phi3 phi4];

end

% QPSK encoding table, Table 17-14.
function phi = qpsk(d_first,d_second)
    table = [ 0; pi/2; pi; 3*pi/2 ];
    phi = table(rebase(d_first*2 + d_second));
end

% get 0-based bit 'i' of 'd'
function b = bit(d,i)
    b = d(:,rebase(i));
end

% rebase 0-based index 'x' to a 1-based index 'y'
function y = rebase(x)
    y = x +1;
end

% [EOF]
