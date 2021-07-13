%Q1 a)
f = 1000; % frequency in Hz.
%generating 1000 samples (even) for 10 cycles- giving about 100 points for each cycle- sufficient accuracy.
%Sampling freq. - 10^5 samples/s
fs = 10^5; %100kHz
t = linspace(0,10*1/f,1000 + 1);
%amplitude 1V sine wave
y = 1*sin(2*pi*f*t);

plot(1000*t,y)
title('Voltage vs time')
xlabel('time -> (ms)')
ylabel('Voltage -> (V)')
grid on
legend('sine wave')
%b)
ydft = fft(y);
N = length(y);
ydft = ydft(1:N/2+1); %[0 or dc  to pi freq.] or [0 to fs/2]
y_enrg = (abs(ydft).^2);
%so if we want the spectral energy in (-pi, pi], we need to twice the samples between 0,pi 
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
legend('sine wave')

%c)
%note- sampling freq 100kHz >> 2*1kHz 
%time window size is fixed 10cycles*1/f = 10ms
%time resolution or sampling period = 1/Fs = 10^-5 s or 0.01 ms
%freq res. is fixed, based on sampling period and no. of time samples
%taken(which is 10*T*sampling freq.)
%frequency resolution- fs/N = 100Hz
%frequency window size - for PSD- 0 to fs/2 or 0 to 50Khz because the signal is real