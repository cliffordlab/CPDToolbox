function [Noisy_Signal] = SNR_Set(Signal, Desired_SNR_dB)
%
%   SNR_Set(x, Desired_SNR_dB) returns the real-valued 
%   input 'Signal' contaminated with normally-distributed, 
%   zero-mean, random noise.  The signal-to-noise ratio   
%   (SNR in dB) of the output 'Noisy_Signal' signal is  
%   controlled by the input 'Desired_SNR_dB' argument measured 
%   in dB.

%   Example:
% 	Npts = 128;                  % Number of time samples
% 	n = 0:Npts-1;                % Time-domain index
% 	Signal = 3*sin(2*pi*3*n/Npts); % Real-valued signal
% 	Desired_SNR_dB = 3;      % Set SNR of output 'Noisy_Signal' to +3 dB
% 	[Noisy_Signal] = SNR_Set(Signal, Desired_SNR_dB);
%
%   Author: Richard Lyons [December 2011]
%******************************************

Npts = length(Signal); % Number of input time samples
Noise = randn(1,Npts); % Generate initial noise; mean zero, variance one

Signal_Power = sum(abs(Signal).*abs(Signal))/Npts;
Noise_Power = sum(abs(Noise).*abs(Noise))/Npts;
	%Initial_SNR = 10*(log10(Signal_Power/Noise_Power));

K = (Signal_Power/Noise_Power)*10^(-Desired_SNR_dB/10);  % Scale factor

New_Noise = sqrt(K)*Noise; % Change Noise level
	%New_Noise_Power = sum(abs(New_Noise).*abs(New_Noise))/Npts
	%New_SNR = 10*(log10(Signal_Power/New_Noise_Power))

Noisy_Signal = Signal + New_Noise;
