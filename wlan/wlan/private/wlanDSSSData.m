function y = wlanDSSSData(PSDU,cfgDSSS)
%wlanDSSSData DSSS processing of the PSDU
%
%   Note: This is an internal undocumented function and its API and/or
%   functionality may change in subsequent releases.
%
%   Y = wlanDSSSData(PSDU,CFGDSSS) generates a DSSS modulated PSDU
%   time-domain waveform for the input PLCP Service Data Unit (PSDU).
%
%   Y is the time-domain PSDU field signal. It is a complex vector of size
%   Ns-by-1, where Ns represents the number of time-domain samples.
%
%   PSDU is the PLCP service data unit input to the PHY. It is a double
%   or int8 typed column vector of length CFGDSSS.PSDULength*8, with each
%   element representing a bit.
%
%   CFGDSSS is the format configuration object of type <a href="matlab:help('wlanNonHTConfig')">wlanNonHTConfig</a>
%   which specifies the parameters for the Non-HT format. Only DSSS
%   modulation type is supported.
%
%   Example: 
%   %  Create a PSDU waveform for 802.11g ERP-CCK 11Mbps operation
%   %  with short preamble:
%   
%      cfgDSSS = wlanNonHTConfig('Modulation','DSSS');
%      cfgDSSS.Preamble = 'Short';
%      cfgDSSS.DataRate = '11Mbps';
%
%      bitsPSDU = randi([0 1],cfgDSSS.PSDULength*8,1);
%
%      waveformPSDU = wlanDSSSData(bitsPSDU,cfgDSSS);
%
%   See also wlanNonHTConfig, wlanDSSSPreamble, wlanDSSSHeader.

%   Copyright 2015 The MathWorks, Inc.

%#codegen
    
    % Only applicable for DSSS modulation of the nonHT format configuration
    dsssValidateConfig(cfgDSSS,mfilename);
    % For DSSS DataRate='1Mbps', Preamble='Short' should be ignored and
    % Preamble='Long' used instead.
    if (strcmpi(cfgDSSS.Preamble,'Short') && strcmpi(cfgDSSS.DataRate,'1Mbps'))
        cfgDSSS.Preamble = 'Long';
    end
    
    % Create information structure
    cfgInfo = dsssInfo(cfgDSSS);
    
    % Validate PSDU input
    validateattributes(PSDU, {'double', 'int8'}, ...
    {'real', 'integer', 'binary', 'column', 'size', [cfgDSSS.PSDULength*8 1]}, ...
    'wlan:wlanDSSSData:InvalidInputPSDU', 'PSDU input');
    
    % Create preamble bits; the preamble bits are scrambled and modulated
    % here to bring the scrambler and modulator to the correct state before
    % processing the header.
    % See dsssInfo(cfgDSSS) for details of SYNC and SFD fields
    preambleBits = [cfgInfo.Sync; cfgInfo.SFD];
            
    % Create header bits; the header bits are scrambled and modulated here
    % to bring the scrambler and modulator to the correct state before
    % processing the PSDU.
    CRC = dsssCRCGenerate([cfgInfo.Signal; cfgInfo.Service; cfgInfo.Length]);
    headerBits = [cfgInfo.Signal; cfgInfo.Service; cfgInfo.Length; CRC];
    
    % Scramble preamble + header + data
    % Clause 17.2.4 PLCP/High Rate PHY data scrambler
    % (same as Clause 16.2.4 PLCP/DSSS PHY data scrambler)
    scrambled = dsssScramble([preambleBits; headerBits; int8(PSDU)], ...
                    cfgInfo.ScramblerInitialization);

    % Modulate preamble + header + data
    % 'P' is the position of the last header symbol in the symbol vector.
    % For DBPSK/DQPSK, the data part will be extracted and spread; for CCK,
    % the last header symbol only shall be used as a reference phase.
    isCCK = any(strcmpi(cfgDSSS.DataRate,{'5.5Mbps','11Mbps'}));
    if (strcmpi(cfgDSSS.DataRate,'1Mbps'))
        
        % Clause 17.2.3.8 Long PLCP data modulation
        % (same as Clause 16.2.5 PLCP data modulation)
        P = length([preambleBits; headerBits]);
        % DBPSK modulation
        pskSymbols = dsssPSKModulate(scrambled,'1Mbps');
        
    else
        
        % Clause 17.2.3.8 Long PLCP data modulation
        % Clause 17.2.3.15 Short PLCP data modulation
        % Repeat appropriate bits to get DBPSK phase
        % transitions through the DQPSK modulator
        if (strcmpi(cfgDSSS.Preamble,'Long'))
            % Repeat preamble and header bits; 'P' is also the right
            % number of modulated symbols to skip when extracting the
            % data part
            P = length([preambleBits; headerBits]);
        else
            % Repeat just preamble bits; for short preamble the header
            % is already DQPSK; 'P' must be amended to skip the
            % modulated header symbols
            P = length(preambleBits);
        end
        % Repeat appropriate bits
        scrambled = [reshape(repmat(scrambled(1:P,1),1,2).',P*2,1); ...
                     scrambled(P+1:end,1)];
        % Amend 'P' for short preamble
        if (strcmpi(cfgDSSS.Preamble,'Short'))
            P = P + length(headerBits)/2;
        end
        % DQPSK modulation
        pskSymbols = dsssPSKModulate(scrambled,'2Mbps');
        
    end
    if (isCCK)
        
        % Extract scrambled PSDU
        scrambledPSDU = scrambled(end-length(PSDU)+1:end,1);
        
        % Extract final header symbol required to manage header to PSDU
        % phase transition and use it to initialise the CCK modulation
        % reference phase
        refPhase = angle(pskSymbols(P));
        
        % Clause 17.4.6.6 Spreading sequences and modulation 
        % for CCK modulation at 5.5 Mb/s and 11 Mb/s
        % CCK modulation
        cckSymbols = dsssCCKModulate(scrambledPSDU,cfgDSSS.DataRate,refPhase);
        
        % Clause 17.4.6.6 Spreading sequences and modulation 
        % for CCK modulation at 5.5 Mb/s and 11 Mb/s
        % CCK spreading
        y = dsssCCKSpread(cckSymbols);
        
    else
        
        % Extract modulated data
        pskSymbols = pskSymbols((P+1):end,1);

        % Spread modulated data
        % Clause 17.4.6.5 Spreading sequence and modulation
        % for 1 Mb/s and 2 Mb/s
        % (same as Clause 16.4.6.4 Spreading sequence)
        % Barker spreading
        y = dsssBarkerSpread(pskSymbols);
        
    end
    
end

% [EOF]