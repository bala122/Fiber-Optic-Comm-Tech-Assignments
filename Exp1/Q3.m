
%PRBS9

% We take 10 PRBS sequences of given order .
% Generating PRBS sequences 
no_seq=10;
y_prbs9= prbs(9,no_seq*9);
y_prbs13= prbs(13,no_seq*13);
y_prbs15= prbs(15,no_seq*15);

% Sampling rate of 100kHz
fs = 10^5;
f = 1000; % frequency of repetition = 1kHz bitrate
% generating time domain pulses from the PRBS sequences
N_symb= fs*1/f;% samples per symb
%replicating value of each sequence over time window.
y_p9 = repmat(transpose(y_prbs9),1,N_symb);
y_p13 = repmat(transpose(y_prbs13),1,N_symb);
y_p15 = repmat(transpose(y_prbs15),1,N_symb);
o1=9;
o2=13;
o3=15;
c1 = [];
c2 = [];
c3 = [];
%time corresponding to 10 sequences
t1 = linspace(0,10*o1*1/f, fs*10*o1*1/f + 1);
t1=t1(1:end-1);
t2 = linspace(0,10*o2*1/f, fs*10*o2*1/f + 1);
t2=t2(1:end-1);
t3 = linspace(0,10*o3*1/f, fs*10*o3*1/f + 1);
t3=t3(1:end-1);
%concatenating and generating the final time domain signal.
for i=1:no_seq*9
    c1 = [ c1 y_p9(i,:)];   
end
for i=1:no_seq*13
    c2 = [ c2 y_p13(i,:)];
end
for i=1:no_seq*15
    c3 = [ c3 y_p15(i,:)];
end


   
plot(1000*t1,c1)
xlabel('time (ms) ->')
ylabel('Voltage (V) ->')
title(' Voltage vs time (PRBS9) ')
grid on
legend('PRBS9')

figure
plot(1000*t2,c2)
xlabel('time (ms) ->')
ylabel('Voltage (V) ->')
title(' Voltage vs time (PRBS13) ')
grid on
legend('PRBS13')

figure
plot(1000*t3,c3)
xlabel('time (ms) ->')
ylabel('Voltage (V) ->')
title(' Voltage vs time (PRBS15)')
grid on
legend('PRBS15')

%b)
ydft1 = fft(c1);
ydft2 = fft(c2);
ydft3 = fft(c3);
N_tot1 = length(c1);
N_tot2 = length(c2);
N_tot3 = length(c3);


ydft1 = ydft1(1:N_tot1/2+1); %[0 or dc  to pi freq.] or [0 to fs/2]
y_enrg1 = (abs(ydft1).^2);
ydft2 = ydft2(1:N_tot2/2+1); %[0 or dc  to pi freq.] or [0 to fs/2]
y_enrg2 = (abs(ydft2).^2);
ydft3 = ydft3(1:N_tot3/2+1); %[0 or dc  to pi freq.] or [0 to fs/2]
y_enrg3 = (abs(ydft3).^2);

%so if we want the spectral energy in (-pi, pi], we need to twice the samples 
%excluding those at 0, pi. (The signal is also real so symmetric conjugate
%FT).
y_enrg1(2:end-1) = 2*y_enrg1(2:end-1);
y_enrg2(2:end-1) = 2*y_enrg2(2:end-1);
y_enrg3(2:end-1) = 2*y_enrg3(2:end-1);
%N time samples, and fs occupies -pi to pi, hence the division to get PSD.
y_psd1 = (y_enrg1/N_tot1)/fs;
y_psd2 = (y_enrg2/N_tot2)/fs;
y_psd3 = (y_enrg3/N_tot3)/fs;

freq1= linspace(0,fs/2,N_tot1/2+1);
freq2= linspace(0,fs/2,N_tot2/2+1);
freq3= linspace(0,fs/2,N_tot3/2+1);

figure
plot(10^-3*freq1,10*log10(y_psd1))
grid on
title('Power Spectral Density vs frequency (PRBS9)')
ylabel('PSD ->( dB/Hz)')
xlabel('Frequency -> (kHz)')
legend('PRBS9')

figure
plot(10^-3*freq2,10*log10(y_psd2))
grid on
title('Power Spectral Density vs frequency (PRBS13)')
ylabel('PSD ->( dB/Hz)')
xlabel('Frequency -> (kHz)')
legend('PRBS13')

figure
plot(10^-3*freq3,10*log10(y_psd3))
grid on
title('Power Spectral Density vs frequency (PRBS15)')
ylabel('PSD ->( dB/Hz)')
xlabel('Frequency -> (kHz)')
legend('PRBS15')

%c)
%sampling rate 100kHz > 2 * 1kHz satisfying nyquist assuming the frequency
%spectrum is bandlimited to the first dip ie, at 1/T(time window of the
%bit) = 1kHz. Ofcourse, we know in reality it isn't so we take a
%sufficiently large sampling rate. Time windows:
%PRBS9:90ms.
%PRBS13:130ms.
%PRBS15:150ms.

%N = 10*order*(1/f)*fs
%freq res = f/(10*order)