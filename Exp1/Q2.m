%Q2.
%a)
f = 1000;
%again we can choose a sampling frequency of approximately 100kHz, because
%ideally, sampling freq >= 2* highest freq. present, in a sq. wave there's
%infinite odd harmonics with decreasing strength, choosing 100Khz we
%satisfy nyquist for upto around 49kHz, which should be sufficient for
%representing the signal.

fs = 10^5; %sampling freq
% time upto 10ms (10 cycles) with appropriate resolution
t = linspace(0,10*1/f, fs* 10*1/f + 1);

%square() has period 2*pi
y = (square(2*pi*f*t, 50)+1)/2;
%plotting in ms
plot(1000*t,y)
xlabel('time (ms) ->')
ylabel('Voltage (V) ->')
title('Square wave Voltage vs time ')
grid on
legend('Square wave (50%)')

%b) 
ydft = fft(y);
N = length(y);
ydft = ydft(1:N/2+1); %[0 or dc  to pi freq.] or [0 to fs/2]
y_enrg = (abs(ydft).^2);
%so if we want the spectral energy in (-pi, pi], we need to twice the samples 
%excluding those at 0, pi. (The signal is also real so symmetric conjugate
%FT).
y_enrg(2:end-1) = 2*y_enrg(2:end-1);
%N time samples, and fs occupies -pi to pi, hence the division to get PSD.
y_psd = (y_enrg/N)/fs;

freq= linspace(0,fs/2,N/2+1);
figure
plot(10^-3*freq,10*log10(y_psd))
grid on
title('Power Spectral Density vs frequency')
ylabel('PSD ->( dB/Hz)')
xlabel('Frequency -> (kHz)')
legend('Square wave (50%)')
%c 
%freq res. = 100Hz (total 1000 samples)