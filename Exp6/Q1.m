%Q1

%PD params
PD_bw = 12*10^9; %Typical BW of PD = 12GHz (above fbaud)
Rd = 0.5; %Typical value of 0.5 A/W
RL = 50; %Typical load resistance
%Defining known constts.
q = 1.6*10^(-19);
k = 1.38064852 *10^(-23);
T = 300;
Fn = 1; %Assuming NF from amplifiers is 0dB 


att = 500; %attenuation of close to ideal LPF

%Generating PRBS


%Input Data Parameters
%Data rate 10GHz
%sps = 4
%No. of symbols 2^12

no_symb = 2^12;
y_prbs9= prbs(9,no_symb);
f_baud = 10*10^9;
t_baud = 1/f_baud;
Nt1 = 4; %sps
fs = Nt1*f_baud;
% Repeating and then concatenating
y_p9_2 = repmat(transpose(y_prbs9),1,Nt1);
prbs_seq=[];
for i=1:no_symb
    prbs_seq = [ prbs_seq y_p9_2(i,:)];   
end

%Power input

%The  " power " = abs(E)^2 calculated from the following value would be
%about 0dBm (peak) which sounds reasonable for a test power given typical
%NEP values of photodetectors .
P_peak = 1*10^(-3);
E_peak = P_peak^0.5;
E_field = prbs_seq*E_peak ;

%PSD of inp E-field calc
N_sig = length(E_field);
freq = linspace(-fs/2, fs/2,N_sig+1);
freq = freq(1:end-1);
[PSD_Ein, freq] = periodogram(E_field,hamming(N_sig),freq,fs);
P_in = (abs(E_field)).^2; %input power


%Noiseless voltage signal generation
Vo1_nn = Rd*RL*((abs(E_field)).^2); 
Vo1_nn = transpose(Vo1_nn); %correcting orientation

%We want signal from -10GHz to +10GHz for avg power
%PSD Calculated after Low pass filtering
Vo1_lpf_nn = lowpass(Vo1_nn,PD_bw,fs,'StopbandAttenuation',att); 
freq = linspace(-fs/2, fs/2,length(Vo1_nn)+1);
freq = freq(1:end-1);
%Calculating PSD
[PSD_vo1_nn, freq] = periodogram(Vo1_lpf_nn,hamming(N_sig),freq,fs);

%Noiseless Voltage PSD output
plot(10^-9*freq,10*log10(PSD_vo1_nn))
xlabel('Frequency (GHz) ->')
ylabel(' PSD dB/(Hz) ->')
title('Power spectral density of noiseless voltage signal')
grid


%Finding Net Signal from photodetector
Vout = PD(E_field, Rd, PD_bw, fs, RL,att); %Noisy volt signal
freq = linspace(-fs/2, fs/2,length(Vout)+1);
freq = freq(1:end-1);

%Finding PSD of the net Voltage Output Signal
[PSD_vnet, freq] = periodogram(Vout,hamming(N_sig),freq,fs);
figure
plot(10^-9*freq,10*log10(PSD_vnet))
xlabel('Frequency (GHz) ->')
ylabel(' PSD dB/(Hz) ->')
title('Power spectral density of output signal(voltage)')
grid



%Calculating Noise using difference of Noiseless, noisy Voltages
V_noise = Vo1_nn - Vout;
%Now that we have the noise, we still restrict it to the PD Bandwitdth
%through the LPF to ensure that our calculations are accurate
V_noise_lpf = lowpass(V_noise,PD_bw,fs,'StopbandAttenuation',att);
%Calculating PSD of noise
PSD_vnoise_lpf = PSD_plot(V_noise_lpf,fs);
figure
plot(10^-9*freq,10*log10(PSD_vnoise_lpf))
xlabel('Frequency (GHz) ->')
ylabel(' PSD (V/m)^2/(Hz) ->')
title('Power spectral density of Noise(voltage)')
grid

df = fs/N_sig;
Pav_in= sum(PSD_Ein) *df
%'avg signal power empirical'
fprintf('power of E-field in via psd :'+string(Pav_in)+'\n');

%Calculating and printing Noise power calculated via PSD

fprintf('noise power via PSD :'+string(sum(PSD_vnoise_lpf)*df)+'\n')

%Theoretical average noise
fprintf('Theoretical avg noise variance :'+string( (2*Rd*mean(P_in)*PD_bw*q*RL^2 +4*k*T*Fn*PD_bw*RL))+'\n')

fprintf('Shot: '+string( (2*Rd*mean(P_in)*PD_bw*q*RL^2) )+'\n')

fprintf('Thermal: '+ string( (4*k*T*Fn*PD_bw*RL))+'\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PSD plot function using FFT
function psd_out = PSD_plot(sig,fs)
%Plotting PSD
ydft = fft(sig);
%Making it two sided from -fs/2 to fs/2
ydft = fftshift(ydft);
N_tot = length(sig);
y_enrg = (abs(ydft).^2);
%Here, we take all the samples to give a 2 sided PSD 
%N time samples, and fs occupies -pi to pi, hence the division to get PSD.
y_psd = (y_enrg/N_tot)/fs;
psd_out = y_psd;


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PHOTODETECTOR TRANSFER FUNCTION
%Ein- column vector
%bw is one sided PD bw.
%NOTE- We assume the Ein is due to transmitted signal + ambient sources

function Vout = PD( Ein, Rd, bw, fs,RL,att)

%Finding power
Pin = (abs(Ein)).^2;

%Defining noiseless o/p current (including dark current)
Id_1 = Rd*Pin;
I_av = mean(Id_1);

%Defining known constts.
q = 1.6*10^(-19);
k = 1.38064852 *10^(-23);
T = 300;
Fn = 1; %Assuming NF from amplifiers is 0dB 

%shot and thermal noise is given by:
%shot: sigma_s^2 =  2qRd*Pin_tot* bw (of PD)
%Thermal: sigma_t^2 = 4kT*Fn*bw/RL 
%We also have to traverse through the  

N = length ( Ein);
Id_net = zeros(N,1);
noise_curr = zeros(N,1);


%Traversing through samples and adding noise accordingly with variance
%having the delta_f of fs/N. N samples would give variance with delta f
%being fs. Now, we add a Low pass filter in the end to restrict it to the
%bandwidth of the PD (bw). Effectively making the delta_f overall to be bw of PD.
for i = 1 :N
sigma_t= sqrt( 2*q*Id_1(i)*fs/N +4*k*T*Fn*fs/(N*RL));
%Mulitplying std dev required to a unit variance white noise
noise_curr(i) = sigma_t*randn;
%Adding white noise of given variance to signal 
Id_net(i) = Id_1(i) + noise_curr(i) ;
end

%Low pass filter with BW = PD one sided BW
Id_net_lpf = lowpass(Id_net,bw,fs,'StopbandAttenuation',att);

%Net output voltage
Vout = Id_net_lpf*RL;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%















