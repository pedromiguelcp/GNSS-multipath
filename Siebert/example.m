clear all;
close all;
clc;
%% load replica
% GPS L1 CA PRN 1 signal and its derivative. The signal is bandlimited to a
% one-sided bandwidth of 10MHz using 20MHz sampling frequency.
load('L1_CA_20MHz.mat', 'signal', 'signal_deriv')


%% shift for corr bank
% Correlators are aligned with samples so we can simply use circshift. In
% this particular case, this corresponds to a correlor spacing of 0.05115
% chip durations
N_corr = 41; % number of correlators
signals = zeros(length(signal), N_corr);
signals_deriv = zeros(length(signal), N_corr);
shiftidx= zeros(N_corr,1);
for iCorr = 1:N_corr
    signals(:, iCorr) = circshift(signal, iCorr - ceil(N_corr/2));
    signals_deriv(:, iCorr) = circshift(signal_deriv, iCorr - ceil(N_corr/2));
    shiftidx(iCorr,1)=iCorr-ceil(N_corr/2);
end

%% Determine autocorrelation functions
phi_ss = signals' * signals;
phi_ss_deriv = signals' * signals_deriv;