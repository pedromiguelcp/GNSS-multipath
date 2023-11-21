function signalShifted = subsampleSignalShifter(signal, tau, FS)
%% Function to perform a temporal shift of a bandlimited signal
% shifting is performed in the frequency domain by a phase manipulation
% enabling subsample shifts
%
%
% Inputs:
% - signal:         the input signal
% - tau_t:          the desired temporal shift in seconds
% - FS:             sampling frequency
%
% Outputs:
% - signalShifted:  signal shifted by tau_t
%
%
% Created by: Christian Siebert, DLR, 03.11.2023

% for compatiability reasons
signal = signal(:);
tau = tau(:);

%% implementation from Pei, Soo-Chang, and Yun-Chiu Lai.
% "Closed form variable fractional delay using FFT with transition band
% trade-off." 2014 IEEE International Symposium on Circuits and Systems
% (ISCAS). IEEE, 2014.
tau_samples = tau*FS;
N = length(signal);
M = length(tau);
tmp = exp(-1i * tau_samples * (2 * pi / N) .* (1:ceil(N/2-1)));
if mod(N, 2) == 0
    h = [ones(M, 1), tmp, cos(tau_samples*pi), conj(flip(tmp, 2))];
else
    h = [ones(M, 1), tmp, conj(flip(tmp, 2))];
end
caSignalFft = fft(signal(:)); % bring signal into frequency domain
caSignalFftShifted = caSignalFft .* transpose(h); % apply frequency ramp

%% bring shifted signal back into time domain
if isreal(signal) % remove imag part when input is real, imag part contains only nummerical errors from matlab
    signalShifted = real(ifft(caSignalFftShifted));
elseif isreal(signal * 1i) % remove real part when input is purely imaginary, real part contains only nummerical errors from matlab
    signalShifted = 1i * imag(ifft(caSignalFftShifted));
else
    signalShifted = ifft(caSignalFftShifted);
end

end
