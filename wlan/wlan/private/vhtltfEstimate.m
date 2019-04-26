function est = vhtltfEstimate(sym,chanBW,nsts,ind)
%vhtltfEstimate Channel estimate using the VHT-LTF
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   EST = vhtltfEstimate(SYM,CHANBW,NSTS,IND) returns channel estimate
%   for each subcarrier specified by the indices IND, using received
%   symbols SYM, channel bandwidth CHANBW, number of space-time streams
%   NUMSTS.

%   Copyright 2015 The MathWorks, Inc.

%#codegen

if (nsts==1)
    % If one space time stream then use LS estimation directly
    ltf = vhtltfSequence(chanBW,nsts);
    est = bsxfun(@rdivide,squeeze(sym(:,1,:)),ltf(ind));
    est = permute(est,[1 3 2]);
else               
    % MIMO channel estimation as per Perahia, Eldad, and Robert Stacey.
    % Next Generation Wireless LANs: 802.11 n and 802.11 ac. Cambridge
    % university press, 2013, page 100, Eq 4.39.
    [ltf,P,nltf] = vhtltfSequence(chanBW,nsts);

    % Verify enough symbols to estimate
    nsym = size(sym,2);
    coder.internal.errorIf(nsym<nltf, ...
        'wlan:wlanChannelEstimate:NotEnoughSymbols',nsts,nltf,nsym);

    Puse = P(1:nsts,1:nltf)'; % Extract and conjugate the P matrix 
    denom = nltf.*ltf(ind);
    nrx = size(sym,3);
    est = complex(zeros(numel(denom),nsts,nrx));
    for i=1:nrx
        rxsym = squeeze(sym(:,(1:nltf),i)); % Symbols on 1 receive antenna
        for j=1:nsts
            est(:,j,i) = rxsym*Puse(:,j)./denom;
        end
    end
end
end