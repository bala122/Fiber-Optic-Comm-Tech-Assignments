%Q_1
%Generating 2^12 prbs symbols using the functions prbs of order 9.
%Directly done because the prbs function repeats the symbols after 2^O -1
%symbols.
no_symb_orig = 2^12;
y_prbs9= prbs(9,no_symb_orig);
f_baud = 10*10^9;
t_baud = 1/f_baud;
%10 Gbaud  is symbol rate => 0.1 ns is symbol duration.
%Since its a square PRBS wave, in discrete time you would need atleast only
%one sample per duration ( sampled at the middle) to represent the whole
%waveform. 
%So we can plot the symbol sequence by taking the samples in the middle of the
%time slot of the OOK 10Gbaud waveform.

%For this we can take a timescale corresponding to only timestamps at the
%middle of the time slot of the OOK waveform.
%So, if we want 10ns, we need 100 symbols since 10ns / 0.1 ns = 100
%timestamps: T-T/2, 2*T-T/2 .... 100*T-T/2.

t = zeros(1,100);
for i=1 :100
    t(i) = i*t_baud -t_baud/2;
end

%Now, we need to truncate y_prbs to 100 symbols.
y_p9_1 = y_prbs9(1:100);

stem(10^9*t,y_p9_1)

xlabel('time (ns) ->')
ylabel('Symbol ->')
title(' Symbols vs time (PRBS9) ')
grid on
legend('PRBS9')


%Q2,3
%However, now we want to generate a continuous waveform, so we
%increase sampling frequency. The sampling rateis 
%Nt = 2, 4, 8  samples per symbol duration
%Hence for 10ns,  total samples =  100 symbols* Nt(samples/symbol)
%sampling rate = fs = Nt (samples/symbol) * fbaud (symbols/s)

%We just plot for 10ns. The original signal as such spans all 2^12 symbols. 
%This original waveform is used for the PSD calculation too.

Nt1 = 2;
no_symb1= 2^12;
fs1 = Nt1*f_baud;
% Repeating and then concatenating
y_p9_2 = repmat(transpose(y_prbs9),1,Nt1);
c1=[];
for i=1:no_symb1
    c1 = [ c1 y_p9_2(i,:)];   
end

% We have to truncate c1 to 10ns duration as well
%10ns corresponds to fs*10*10^-9 samples

c1_copy = c1(1:fs1*10*10^-9);


%plotting for 10ns.
t1 = linspace(0,10*10^-9, fs1*10*10^-9 + 1);
t1=t1(1:end-1);
figure
plot(10^9*t1,c1_copy)
xlabel('time (ns) ->')
ylabel('Voltage (V) ->')
title(' Voltage vs time (PRBS9, Nt=2) ')
grid on
legend('PRBS9(Nt=2)')

%Plotting PSD- This is using all the symbols there.
ydft1 = fft(c1);
N_tot1 = length(c1);
ydft1 = ydft1(1:N_tot1/2+1); %[0 or dc  to pi freq.] or [0 to fs/2]
y_enrg1 = (abs(ydft1).^2);
%so if we want the spectral energy in (-pi, pi], we need to twice the samples 
%excluding those at 0, pi. (The signal is also real so symmetric conjugate
%FT).
y_enrg1(2:end-1) = 2*y_enrg1(2:end-1);
%N time samples, and fs occupies -pi to pi, hence the division to get PSD.
y_psd1 = (y_enrg1/N_tot1)/fs1;
freq1= linspace(0,fs1/2,N_tot1/2+1);

figure
plot(10^-9*freq1,10*log10(y_psd1))
grid on
title('Power Spectral Density vs frequency (PRBS9, Nt =2)')
ylabel('PSD ->( dB/Hz)')
xlabel('Frequency -> (GHz)')
legend('PRBS9 (Nt=2)')


%Nt = 4


Nt2 = 4;
no_symb2= 2^12;
fs2 = Nt2*f_baud;
% Repeating and then concatenating
y_p9_3 = repmat(transpose(y_prbs9),1,Nt2);
c2=[];
for i=1:no_symb2
    c2 = [ c2 y_p9_3(i,:)];   
end

% We have to truncate c2 to 10ns duration as well
%10ns corresponds to fs*10*10^-9 samples

c2_copy = c2(1:fs2*10*10^-9);

%Plotting for 10ns.
t2 = linspace(0,10*10^-9, fs2*10*10^-9 + 1);
t2=t2(1:end-1);
figure
plot(10^9*t2,c2_copy)
xlabel('time (ns) ->')
ylabel('Voltage (V) ->')
title(' Voltage vs time (PRBS9, Nt=4) ')
grid on
legend('PRBS9(Nt=4)')


%Plotting PSD
ydft2 = fft(c2);
N_tot2 = length(c2);
ydft2 = ydft2(1:N_tot2/2+1); %[0 or dc  to pi freq.] or [0 to fs/2]
y_enrg2 = (abs(ydft2).^2);
%so if we want the spectral energy in (-pi, pi], we need to twice the samples 
%excluding those at 0, pi. (The signal is also real so symmetric conjugate
%FT).
y_enrg2(2:end-1) = 2*y_enrg2(2:end-1);
%N time samples, and fs occupies -pi to pi, hence the division to get PSD.
y_psd2 = (y_enrg2/N_tot2)/fs2;
freq2= linspace(0,fs2/2,N_tot2/2+1);

figure
plot(10^-9*freq2,10*log10(y_psd2))
grid on
title('Power Spectral Density vs frequency (PRBS9, Nt =4)')
ylabel('PSD ->( dB/Hz)')
xlabel('Frequency -> (GHz)')
legend('PRBS9 (Nt=4)')




%Nt = 8


Nt3 = 8;
no_symb3= 2^12;
fs3 = Nt3*f_baud;
% Repeating and then concatenating
y_p9_4 = repmat(transpose(y_prbs9),1,Nt3);
c3=[];
for i=1:no_symb3
    c3 = [ c3 y_p9_4(i,:)];   
end

% We have to truncate c3 to 10ns duration as well
%10ns corresponds to fs*10*10^-9 samples

c3_copy = c3(1:fs3*10*10^-9);

%plotting for 10ns.
t3 = linspace(0,10*10^-9, fs3*10*10^-9 + 1);
t3=t3(1:end-1);
figure
plot(10^9*t3,c3_copy)
xlabel('time (ns) ->')
ylabel('Voltage (V) ->')
title(' Voltage vs time (PRBS9, Nt=8) ')
grid on
legend('PRBS9(Nt=8)')

%Plotting PSD
ydft3 = fft(c3);
N_tot3 = length(c3);
ydft3 = ydft3(1:N_tot3/2+1); %[0 or dc  to pi freq.] or [0 to fs/2]
y_enrg3 = (abs(ydft3).^2);
%so if we want the spectral energy in (-pi, pi], we need to twice the samples 
%excluding those at 0, pi. (The signal is also real so symmetric conjugate
%FT).
y_enrg3(2:end-1) = 2*y_enrg3(2:end-1);
%N time samples, and fs occupies -pi to pi, hence the division to get PSD.
y_psd3 = (y_enrg3/N_tot3)/fs3;
freq3= linspace(0,fs3/2,N_tot3/2+1);

figure
plot(10^-9*freq3,10*log10(y_psd3))
grid on
title('Power Spectral Density vs frequency (PRBS9, Nt =8)')
ylabel('PSD ->( dB/Hz)')
xlabel('Frequency -> (GHz)')
legend('PRBS9 (Nt=8)')


%Q4,5
%No noise
%scatterplot parameters 'o', 1,0 are just the marker params/ offset etc.
scatterplot(y_prbs9,1,0,'o')
grid on
xlabel('I-phase')
ylabel('Q-phase')
title('Constellation diagram for OOK(no noise)')
legend('OOK (no noise)')

%Now, adding noise using the awgn function after measuring input power.
%Note- Here noise is added to each symbol which is plotted on the
%constellation diagram for it's SNR.
%5dB
y_5dB = awgn(y_prbs9, 5, 'measured');
scatterplot(y_5dB,1,0,'o')
grid on
xlabel('I-phase')
ylabel('Q-phase')
title('Constellation diagram for OOK(SNR 5dB)')
legend('OOK (SNR 5dB)')

%10dB
y_10dB = awgn(y_prbs9(1:2000), 10, 'measured');
scatterplot(y_10dB,1,0,'o')
grid on
xlabel('I-phase')
ylabel('Q-phase')
title('Constellation diagram for OOK(SNR 10dB)')
legend('OOK (SNR 10dB)')

%15dB
y_15dB = awgn(y_prbs9, 15, 'measured');
scatterplot(y_15dB,1,0,'o')
grid on
xlabel('I-phase')
ylabel('Q-phase')
title('Constellation diagram for OOK(SNR 15dB)')
legend('OOK (SNR 15dB)')



%Q6 Eye diagram

%Here, we take the Nt = 8 samples per symbols case for the electrical
%waveform because it's most accurate PRBS representation of the three cases.
%For the eye diagram, we need an electrical waveform, hence the the noise added
%will be to the time domain waveform with a given Nt( samples/ symbol) 
% after measuring input power. 
 
sig5=zeros(length(c3),1);
sig5 = transpose(c3);

out_5dB = awgn(sig5, 5, 'measured');
out_10dB = awgn(sig5,10,'measured');
out_15dB = awgn(sig5,15, 'measured');




%5dB
%input signal for this is out_5dB
%the second parameter 'n' denotes it plots n+1 samples 1st and (n+1)th
%sample being same. So it takes the appropriate number or samples, symbol
%periods.
%3rd parameter denotes total time span (here 0.2ns for 2 symbols ), 
%the time axis is scaled b/w -0.1 ns to 0.1 ns
%last parameter denotes offset in terms of samples.
%Here, Nt = 8, but we take more samples and apply an offset 
%to get the eye diagram in the position.
eyediagram(out_5dB,16,0.2,3)
grid on 
title('Eye diagram OOK(SNR 5dB)')
legend('OOK(SNR 5dB)')
xlabel('Time (ns)')

%Similarly for-
%10dB
eyediagram(out_10dB,16,0.2,3)
grid on 
title('Eye diagram OOK(SNR 10dB)')
legend('OOK(SNR 10dB)')
xlabel('Time (ns)')

%15dB
eyediagram(out_15dB,16,0.2,3)
grid on 
title('Eye diagram OOK(SNR 15dB)')
legend('OOK(SNR 15dB)')
xlabel('Time (ns)')

