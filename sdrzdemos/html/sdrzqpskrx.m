%% QPSK Receiver with Analog Devices FMCOMMS1
%
% This example shows how to use the Xilinx(R) Zynq-Based Radio Support
% Package and Communications System Toolbox(TM) software to implement a
% QPSK receiver in Simulink(R). The receiver addresses practical issues in
% wireless communications, such as carrier frequency and phase offset,
% timing offset and frame synchronization. This model receives the signal
% sent by the <matlab:showdemo('sdrzqpsktx') QPSK Transmitter with SDR
% Hardware> model. The receiver demodulates the received symbols and
% outputs a simple message to the MATLAB(R) command line.
% 
% Refer to the <matlab:sdrzdoc('sdrzspsetup') Getting
% Started> documentation for details on configuring your host computer to
% work with the Support Package for Xilinx(R) Zynq-Based Radio.

% Copyright 2014-2015 The MathWorks, Inc.

%% Introduction
%
% This example receives a QPSK signal over the air using SDR hardware. It
% has two main objectives:
%
% * Implement a prototype QPSK-based receiver in Simulink using the
% Xilinx(R) Zynq-Based Radio Support Package.
%
% * Illustrate the use of key Communications System Toolbox(TM) System
% objects for QPSK system design.
%
% In this example, the <matlab:sdrzdoc('commsinternalsdrrxzc706fmc1sl')
% ZC706 and Analog Devices(TM) FMCOMMS1 Receiver> block receives a signal impaired by the over-the-air
% transmission and outputs complex baseband signals that are processed in
% Simulink. This example provides a sample design of a practical digital
% receiver that can cope with wireless channel impairments. The receiver
% includes:
% 
% * FFT-based coarse frequency compensation
% * PLL-based fine frequency compensation
% * Timing recovery with fixed-rate resampling
% * Bit stuffing/skipping
% * Frame synchronization
% * Phase ambiguity correction


%% Setup
% 
% Before running the example, ensure you have performed the following
% steps:
%
% 1. Configure your host computer to work with the Support Package for
% Xilinx(R) Zynq-Based Radio. See <matlab:sdrzdoc('sdrzspsetup')
% Getting Started> for help.
%
% * Some additional steps may be required if you want to run two radios
% from a single host computer. See
% <matlab:sdrzdoc('sdrz_tworadios') Setup For Two Radios-One Host>
% for help.
%
% 2. Ensure that a suitable signal is available for reception. This example
% is designed to work in conjunction with any of the following possible
% signal sources:
%
% * The <matlab:showdemo('sdrzQPSKTransmitter') QPSK Transmitter with Analog Devices FMCOMMS1> MATLAB example.
% * The <matlab:showdemo('sdrzqpsktx') QPSK Transmitter with Analog Devices FMCOMMS1> Simulink example.
% * The <matlab:showdemo('sdrzQPSKTxFPGA') Targeting HDL-Optimized QPSK Transmitter with Analog Devices FMCOMMS1> Simulink example.


%% Running the Example
%
% Start the transmitter before running the model. When you run the
% model, the received messages are decoded and printed out in the
% MATLAB command window while the simulation is running. If the received
% signal is decoded correctly, you should see 'Hello world 0##' messages in
% the MATLAB command line similar to those shown below.
% 
%  Hello world 031
%  Hello world 032
%  Hello world 033
%  Hello world 034
%  Hello world 035
%  Hello world 036
%  Hello world 037
%  Hello world 038
%  Hello world 039
%  Hello world 040

%%
% You can configure various receiver algorithm parameters using the *Model
% Parameters* block.


%% Receiver Design: System Architecture
%
% The top-level structure of the model is shown in the following figure.
%
modelname = 'sdrzqpskrx';
open_system(modelname);
set_param(modelname, 'SimulationCommand', 'update')
%%
% The detailed structures of the *QPSK Receiver* subsystem are illustrated
% in the following figure.
%%
%
open_system([modelname '/QPSK Receiver']);
%%
% The components are further described in the following sections:
%
% * *Automatic Gain Control* - Applies a variable gain to try keep the
% signal amplitude at _1/Upsampling Factor_
% * *Raised Cosine Receive Filter* - Uses a roll-off factor of 0.5, and
% downsamples the input signal by two
% * *Coarse Frequency Compensation* - Estimates an approximate frequency
% offset of the received signal and corrects it
% * *Fine Frequency Compensation* - Compensates for the residual frequency
% offset and the phase offset
% * *Timing Recovery* - Resamples the input signal according to a recovered
% timing strobe so that symbol decisions are made at the optimum sampling
% instants
% * *Data Decoding* - Aligns the frame boundaries, resolves the phase
% ambiguity caused by the Fine Frequency Compensation subsystem,
% demodulates the signal, and decodes the text message

%% 
% *Automatic Gain Control AGC*
%
% The phase error detector gain $K_p$ of the phase and timing error
% detectors is proportional to the received signal amplitude and the
% average symbol energy. To ensure an optimum loop design, the signal
% amplitude at the inputs of the carrier recovery and timing recovery loops
% must be stable. The AGC sets its output power to _1/Upsampling Factor_
% (0.25),  so that the equivalent gains of the phase and timing error
% detectors remain constant over time. The AGC is placed before the *Raised
% Cosine Receive Filter* so that the signal amplitude can be measured with
% an oversampling factor of four, thus improving the accuracy of the
% estimate. The *AGC* subsystem  updates its compensation gain after each
% block of ten QPSK symbols in order to smooth out variations in signal
% amplitude. You can refer to Chapter 7.2.2 and Chapter 8.4.1 of [ <#19 1>
% ] for details on how to design the phase detector gain $K_p$.

%%
% *Raised Cosine Receive Filter*
%
% The *Raised Cosine Receive Filter* downsamples the input signal by a
% factor of two, with a roll-off factor of 0.5. It provides matched
% filtering for the transmitted waveform.

%%
% *Coarse Frequency Compensation*
%
% The *Coarse Frequency Compensation* subsystem corrects the input signal
% with a rough estimate of the frequency offset. The following diagram
% shows the *Find Frequency Offset* subsystem in the *Coarse Frequency
% Compensation* subsystem. This subsystem uses a baseband QPSK signal with
% a designated phase index $n$, frequency offset $\Delta f$ and phase
% offset $\Delta \phi$ expressed as $e^{(j(n\pi/2+\Delta ft+\Delta
% \phi))}$, $n=0,1,2,3$. First, the subsystem raises the input signal to
% the power of four to obtain $e^{(j(4\Delta ft+4\Delta \phi))}$, which is
% not a function of the QPSK modulation. Then it performs an FFT on the
% modulation-independent signal to estimate the tone at four times the
% frequency offset. After dividing the estimate by four, the
% *Phase/Frequency Offset* library block corrects the frequency offset.
% There is usually a residual frequency offset even after the coarse
% frequency compensation, which would cause a slow rotation of the
% constellation.  The *Fine Frequency Compensation* subsystem compensates
% for this residual frequency.
%
close_system([modelname '/QPSK Receiver']);
open_system([modelname '/QPSK Receiver/Coarse Frequency Compensation/Find Frequency Offset'],'force');

%% 
% *Fine Frequency Compensation*
%
% The Fine Frequency Compensation subsystem implements a phase-locked loop
% (PLL), described in Chapter 7 of [ <#19 1> ], to track the residual
% frequency offset and the phase offset in the input signal, as shown in
% the following figure. The PLL uses a *Direct Digital Synthesizer (DDS)*
% to generate the compensating phase that offsets the residual frequency
% and phase offsets. The phase offset estimate from *DDS* is the integral
% of the phase error output of the *Loop Filter*.
%
% A maximum likelihood *Phase Error Detector (PED)* , described in Chapter
% 7.2.2 of [ <#19 1> ], generates the phase error. A tunable
% proportional-plus-integral *Loop Filter* , described in Appendix C.2 of [
% <#16 1> ] filters the error signal and then feeds it into the *DDS*. The
% _Loop Bandwidth_ (normalized by the sample rate) and the _Loop Damping
% Factor_ are tunable for the *Loop Filter*. The default normalized loop
% bandwidth is set to 0.06 and the default damping factor is set to 2.5
% (over damping) so that the PLL quickly locks to the intended phase while
% introducing minimal phase noise.
%
close_system([modelname '/QPSK Receiver/Coarse Frequency Compensation/Find Frequency Offset']);
open_system([modelname '/QPSK Receiver/Fine Frequency Compensation']);

%%
% *Timing Recovery*
%
% The *Timing Recovery* subsystem implements a PLL, described in Chapter 8
% of [ <#19 1> ], to correct the timing error in the received signal. The
% input of the *Timing Recovery* subsystem is oversampled by two. On
% average the *Timing Recovery* subsystem generates one output sample for
% every two input samples. The *NCO Control* subsystem implements a
% decrementing modulo-1 counter, described in Chapter 8.4.3 of [ <#19 1> ],
% to generate the control signal for the *Modified Buffer*, that selects
% the interpolants of the *Interpolation Filter*. This control signal also
% enables the *Timing Error Detector (TED)*, so that it calculates the
% timing errors at the correct timing instants. The *NCO Control* subsystem
% updates the timing difference for the *Interpolation Filter* , generating
% interpolants at optimum sampling instants. The *Interpolation Filter* is
% a Farrow parabolic filter with $\alpha=0.5$ as described in Chapter 8.4.2
% of [ <#19 1> ]. The filter uses an $\alpha$ of 0.5 so that all the filter
% coefficients become only 1, -1/2 and 3/2, which significantly simplifies
% the interpolator structure. Based on the interpolants, timing errors are
% generated by a zero-crossing *Timing Error Detector*, described in
% Chapter 8.4.1 of [ <#19 1> ], filtered by a tunable
% proportional-plus-integral *Loop Filter*, described in Appendix C.2 of [
% <#19 1> ], and fed into the *NCO Control* for a timing difference update.
% The _Loop Bandwidth_ (normalized by the sample rate) and _Loop Damping
% Factor_ are tunable for the *Loop Filter*. The default normalized loop
% bandwidth is set to 0.01 and the default damping factor is set to unity
% (critical damping) so that the PLL quickly locks to the correct timing
% while introducing little phase noise.
%
% When the timing error (delay) reaches symbol boundaries, there will be
% one extra or missing interpolant in the output.  The TED implements bit
% stuffing/skipping to handle the extra/missing interpolants. You can refer
% to Chapter 8.4.4 of [ <#19 1> ] for details of bit stuffing/skipping.
% 
close_system([modelname '/QPSK Receiver/Fine Frequency Compensation']);
open_system([modelname '/QPSK Receiver/Timing Recovery/Timing Recovery PLL']);
%%
% The timing recovery loop normally generates 100 QPSK symbols per frame,
% one output symbol for every two input samples. It also outputs a timing
% strobe that runs at the input sample rate. Under normal circumstances,
% the strobe value is simply a sequence of alternating ones and zeroes.
% However, this occurs only when the relative delay between Tx and Rx
% contains some fractional part of one symbol period and the integer part
% of the delay (in symbols) remains constant. If the integer part of the
% relative delay changes, the strobe value can have two consecutive zeroes
% or two consecutive ones. In that case, the timing recovery loop generates
% 99 or 101 QPSK output symbols per frame. However, the downstream
% processing must use a frame size of 100 symbols, which is ensured by the
% *Modified Buffer* subsystem.
% 
% The *Modified Buffer* subsystem uses the strobe to fill up a delay line
% with properly sampled QPSK symbols. As each QPSK symbol is added to the
% delay line, a counter increments the number of symbols in the line. At
% each sampling instant, the delay line outputs a frame of size 100 to the
% *Data Decoding* subsystem. However, the *Data Decoding* subsystem runs on
% its received data only when its enable signal goes high. This occurs when
% both the counter value reaches 100 and the strobe is high, i.e. each time
% exactly 100 valid QPSK symbols are present at the *Modified Buffer*.

%% 
% *Data Decoding*
%
% The *Data Decoding* subsystem performs frame synchronization, phase
% ambiguity resolution, demodulation and text message decoding. The
% subsystem uses a QPSK-modulated Barker code, generated by the *Bits
% Generation* subsystem, to correlate against the received QPSK symbols and
% achieve frame synchronization. The *Compute Delay* subsystem correlates
% the data input with the QPSK modulated Barker code, and uses the index of
% the peak amplitude to find the delay.
%
% The carrier phase PLL of the *Fine Frequency Compensation* subsystem may
% lock to the unmodulated carrier with a phase shift of 0, 90, 180, or 270
% degrees, which can cause a phase ambiguity. For details of phase
% ambiguity and its resolution, refer to Chapter 7.2.2 and 7.7 in [
% <#16 1> ]. The *Phase Offset Estimator* subsystem determines this phase
% shift. The *Phase Ambiguity Correction & Demodulation* subsystem rotates
% the input signal by the estimated phase offset and demodulates the
% corrected data. The payload bits are descrambled, and the first 105
% payload bits are extracted and stored in a workspace variable. All the
% stored bits are converted to characters and printed out at the MATLAB
% command window while the simulation is running.


%% Things to Try
%
% The example allows you to experiment with multiple system capabilities to
% examine their effect on bit error rate performance. You can make changes
% to the receiver algorithms behaviors using the *Model Parameters* block.
% 
% * You can tune the _FFT Size_ and _Number of Spectrum Averages_ for the
% *Coarse Frequency Compensation* subsystem to see the effect of the
% estimation accuracy and the tolerance to a high noise level. The
% resolution of the estimate is the frequency spacing between two adjacent
% FFT points. There is a speed versus accuracy tradeoff when choosing the
% value of _FFT Size_. To get a more accurate frequency estimate usually
% requires a larger _FFT Size_. However, a larger _FFT Size_ also incurs a
% higher computational burden. If the resolution of the *Coarse Frequency
% Compensation* subsystem is low, then the *Fine Frequency Compensation*
% subsystem must have a wider frequency tracking range.
% * Due to the existence of noise and zero padding of the input, the FFT
% output might have some outliers in the estimation results. To ease the
% effect of these bad estimates, you can adjust the _Number of Spectrum
% Averages_ to average the FFT result across multiple frames. The larger
% _Number of Spectrum Averages_ improves the robustness of the coarse
% frequency estimation, but this also incurs a greater computational
% burden. Also, the fourth-power operation can correctly estimate an offset
% only if the offset satisfies the following inequality:
% 
% $4*\Delta f_{\max} \le f_s/2$, or
% 
% $4*\Delta f_{\max} \le 2*R_{sym}/2$, or
% 
% $\Delta f_{\max} \le R_{sym}/4$.
% 
% * The FFT-based *Coarse Frequency Compensation* subsystem was
% designed for a scenario with a static frequency offset. In practice, the
% frequency offset might vary over time. This model can still track a
% time-varying frequency drift by the *Coarse Frequency Compensation*
% subsystem. However, the coarse frequency estimates take on discrete
% values, separated by the frequency resolution of the subsystem. You might
% observe jumps between frequency estimates. You can also implement coarse
% frequency compensation with a filter to get a smoother estimation output.
% * You can adjust the PLL design parameters such as _Loop Bandwidth_ and
% _Damping Factor_ in both *Fine Frequency Compensation* and *Timing
% Recovery* subsystems to see their effect on pull-in range, convergence
% time and the estimation accuracy. With a large _Loop Bandwidth_ and
% _Damping Factor_, the PLL can acquire over a greater frequency offset
% range. However a large _Loop Bandwidth_ allows more noise, which leads to
% a large mean squared error in the phase estimation. "Underdamped systems
% (with Damping Factor less than one) have a fast settling time, but
% exhibit overshoot and oscillation; overdamped systems (with Damping
% Factor greater than one) have a slow settling time but no oscillations."
% [ <#19 1> ]. For more detail on the design of these PLL parameters, you
% can refer to Appendix C in [ <#19 1> ].
% * The *Timing Recovery* subsystem relies on a stable constellation which
% does not rotate over time. This requires an accurate frequency offset
% compensation. In this model, if the actual frequency offset exceeds the
% maximum frequency offset that can be tracked by the current coarse
% compensation subsystem, you can increase its tracking range by increasing
% the oversampling factor. Another way to adjust the tracking range is to
% implement a rotationally-invariant timing error detector (e.g., Gardner
% timing error detector described in Chapter 8.4.1 of [ <#19 1> ]) first
% and correct the rotation afterwards.
close_system([modelname '/QPSK Receiver/Timing Recovery/Timing Recovery PLL']);
close_system(modelname, 0);


%% Alternative Implementations
%
% This example describes the Simulink implementation of a QPSK receiver
% with SDR Hardware. You can also view a MATLAB implementation of this
% example in <matlab:showdemo('sdrzQPSKReceiver') QPSK Receiver with SDR
% Hardware using MATLAB>.
%
% You can also explore a non-hardware QPSK transmitter and receiver example
% that models a general wireless communication system using an AWGN channel
% and simulated channel impairments with
% <matlab:showdemo('commQPSKTransmitterReceiver')
% commQPSKTransmitterReceiver>.


%% Troubleshooting the Example
%
% If you fail to successfully receive any 'Hello world' messages, try the
% troubleshooting steps below:
% 
% * If the example runs slower than real time, you can try using
% <matlab:sdrzdoc('sdrz_burstmode') burst mode>.
% * The ability to decode the received signal depends on the received
% signal strength. If the message is not properly decoded by the receiver
% system, you can vary the gain applied to the received signal by changing
% the _Gain_ parameter in the SDR receiver block.
% * A large relative frequency offset between the transmit and receive
% radios can prevent the receiver from properly decoding the message.  If
% that happens, you can determine the offset by sending a tone at a known
% frequency from the transmitter to the receiver, and then measuring the
% offset between the transmitted and received frequency. This value can
% then be used to compensate the center frequency of the receiver block.
% See the <matlab:showdemo('sdrzfreqcalib') Frequency Offset Calibration
% with Analog Devices FMCOMMS1> example.
% 
% If you still fail to receive any messages, see
% <matlab:sdrzdoc('sdrz_troubleshoot') Xilinx Zynq-Based Radio Processing
% Errors and Fixes>.


%% List of Example Helper Files
%
% This example uses the following helper files:
%
% * <matlab:edit('sdrzqpskrx_init') sdrzqpskrx_init.m>: returns a structure
% of variables used to control constant parameters in the model.
% * <matlab:edit('sdrzqpskrx_mask') sdrzqpskrx_mask.m>: controls the *Model
% Parameters* mask and updates the model depending on the parameters
% supplied in the mask.


%% References
%
% 1. Michael Rice, "Digital Communications - A Discrete-Time Approach",
% Prentice Hall, April 2008.

