classdef QPSKBitsGeneratorR < matlab.System
    properties (Nontunable)
        MessageLength = 105;
        BernoulliLength = 69;
        ScramblerBase = 2;
        ScramblerPolynomial = [1 1 1 0 1];
        ScramblerInitialConditions = [0 0 0 0];

        % LDPC
        LdpcNewH = zeros(148,296);
        LdpcU = zeros(148,148);
        LdpcL = zeros(148,148);
    end

    properties (Access=private)
        pHeader
        pScrambler
        pMsgStrSet
        pCount
    end

    methods 
        function obj = QPSKBitsGeneratorR(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end

    methods (Access=protected)
        function setupImpl(obj,~)
            bbc = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];
            bbc = [bbc bbc];
            ubc = ((bbc + 1) /2)';
            temp = (repmat(ubc,1,2))';
            obj.pHeader = temp(:);
            obj.pCount = 0;
            obj.pScrambler = comm.pScrambler(obj.ScramblerBase,...
                                obj.ScramblerPolynomial,...
                                obj.ScramblerInitialConditions);
            obj.pMsgStrSet = ['Hello world 0000';...
            'Hello world 0001';...
            'Hello world 0002';...
            'Hello world 0003';...
            'Hello world 0004';...
            'Hello world 0005';...
            'Hello world 0006';...
            'Hello world 0007';...
            'Hello world 0008';...
            'Hello world 0009';...
            'Hello world 0010';...
            'Hello world 0011';...
            'Hello world 0012';...
            'Hello world 0013';...
            'Hello world 0014';...
            'Hello world 0015';...
            'Hello world 0016';...
            'Hello world 0017';...
            'Hello world 0018';...
            'Hello world 0019';...
            'Hello world 0020';...
            'Hello world 0021';...
            'Hello world 0022';...
            'Hello world 0023';...
            'Hello world 0024';...
            'Hello world 0025';...
            'Hello world 0026';...
            'Hello world 0027';...
            'Hello world 0028';...
            'Hello world 0029';...
            'Hello world 0030';...
            'Hello world 0031';...
            'Hello world 0032';...
            'Hello world 0033';...
            'Hello world 0034';...
            'Hello world 0035';...
            'Hello world 0036';...
            'Hello world 0037';...
            'Hello world 0038';...
            'Hello world 0039';...
            'Hello world 0040';...
            'Hello world 0041';...
            'Hello world 0042';...
            'Hello world 0043';...
            'Hello world 0044';...
            'Hello world 0045';...
            'Hello world 0046';...
            'Hello world 0047';...
            'Hello world 0048';...
            'Hello world 0049';...
            'Hello world 0050';...
            'Hello world 0051';...
            'Hello world 0052';...
            'Hello world 0053';...
            'Hello world 0054';...
            'Hello world 0055';...
            'Hello world 0056';...
            'Hello world 0057';...
            'Hello world 0058';...
            'Hello world 0059';...
            'Hello world 0060';...
            'Hello world 0061';...
            'Hello world 0062';...
            'Hello world 0063';...
            'Hello world 0064';...
            'Hello world 0065';...
            'Hello world 0066';...
            'Hello world 0067';...
            'Hello world 0068';...
            'Hello world 0069';...
            'Hello world 0070';...
            'Hello world 0071';...
            'Hello world 0072';...
            'Hello world 0073';...
            'Hello world 0074';...
            'Hello world 0075';...
            'Hello world 0076';...
            'Hello world 0077';...
            'Hello world 0078';...
            'Hello world 0079';...
            'Hello world 0080';...
            'Hello world 0081';...
            'Hello world 0082';...
            'Hello world 0083';...
            'Hello world 0084';...
            'Hello world 0085';...
            'Hello world 0086';...
            'Hello world 0087';...
            'Hello world 0088';...
            'Hello world 0089';...
            'Hello world 0090';...
            'Hello world 0091';...
            'Hello world 0092';...
            'Hello world 0093';...
            'Hello world 0094';...
            'Hello world 0095';...
            'Hello world 0096';...
            'Hello world 0097';...
            'Hello world 0098';...
            'Hello world 0099']; 
        end

        function [y,msg] = stepImpl(obj)
            cycle = mod(obj.pCount,100);
            msgStr = obj.pMsgStrSet(cycle+1,:);
            msgBin = de2bi(int8(msgStr),7,'left-msb');
            msg = reshape(double(msgBin).', obj.MessageLength,1);
            data = [msg; randi([0 1], obj.BernoulliLength/2,1)];

            % Scramble the data
            scrambledData = step(obj.pScrambler,data);

            lpdccode = mod(obj.LdpcU\(obj.LdpcL\mod(obj.LdpcNewH(:,148+1:end)*scrambledData,2)),2);
            CodescrambledData = [lpdccode;scrambledData];

            y = [obj.pHeader; CodescrambledData];

            obj.pCount = obj.pCount+1;
        end

        function resetImpl(obj)
            obj.pCount = 0;
            reset(obj.pScrambler);
        end

        function releaseImpl(obj)
            release(obj.pScrambler);
        end

        function N = getNumInputsImpl(~)
            N = 0;
        end

        function N = getNumOutputsImpl(~)
            N = 2;
        end
    end
end

