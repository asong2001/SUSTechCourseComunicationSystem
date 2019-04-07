% QPSK tranmitter and receiver 
% intro LDPC for comparison 
prmQPSKTxRx = commqpsktxrx_init %#ok<*NOPTS> % QPSK system parameters

useScopes = true;   % true if scopes are to be used
printReceivedData = false; %true if the received data is to be printed
compileIt = false;  % true if code is to be compiled
useCodegen = false; % true to run the generated mex file

%% Execution and Results
%
% 1) Constellation diagram of the *Raised Cosine Receive Filter* output.
%
% 2) Constellation diagram of the *Symbol Synchronizer* output.
%
% 3) Constellation diagram of the *Fine Frequency Compensation* output.
%
% 4) Power spectrum of the *Raised Cosine Receive Filter* output.

if compileIt
    codegen -report runQPSKSystemUnderTest.m -args {coder.Constant(prmQPSKTxRx),coder.Constant(useScopes),coder.Constant(printReceivedData)} %#ok
end
if useCodegen
    BER = runQPSKSystemUnderTest_mex(prmQPSKTxRx, useScopes, printReceivedData);
else
    BER = runQPSKSystemUnderTest(prmQPSKTxRx, useScopes, printReceivedData);
end
fprintf('Error rate = %f.\n',BER(1));
fprintf('Number of detected errors = %d.\n',BER(2));
fprintf('Total number of compared samples = %d.\n',BER(3));

%% References
% 1. Rice, Michael. _Digital Communications - A Discrete-Time
% Approach_. 1st ed. New York, NY: Prentice Hall, 2008.

displayEndOfDemoMessage(mfilename)
