%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NOTE- The noise addition is a random process and since the plots plotted
%are including dispersive effects the differences between them are minor
%and may sometimes not be as expected, hence running the code again would
%give a better result.

%Q4
%Using APD of gain 1000 case
%Generating OOK mod. at 20Gbaud 
sps = 4;
att =60; %attenuation in LPF stopband

%MZM data
Vpi = 3;
V_off = 0.5;
IL = 6;
ER = 20;
Pin = 0; %in dBm
%Converting to W
Pin = 10^(Pin/10);
Pin = Pin *10^(-3);

%data sequence 
% PRBS9 sps=4 Bitrate 20Gbps
%Hence bit slot time is 0.05ns
%Sampling frequency is fixed as sps*baudrate = 80GHz.
fbaud = 20*10^9;
%Generate the Time domain modulating waveform for 2^12 bits =( 204.8ns)
no_symb= 2^12;
t_duration = no_symb/fbaud;


%Making the time domain input OOK voltage signal by
%repeating the symbol and making it look like a square wave.
y_prbs9= prbs(9,no_symb);
y_p9 = repmat(transpose(y_prbs9),1,sps);
V_sig=[];
for i=1:no_symb
    V_sig = [ V_sig y_p9(i,:)];   
end

%Choosing a swing of Vpi (peak-peak)
V_sig = Vpi*(V_sig- 0.5);
%Biasing @ 0.5 Vpi 
V_sig = V_off+0.5*Vpi+V_sig;
V_sig = transpose(V_sig); % Ein, V_sig are column vectors



%Laser data
%LW 100MHz
lw= 100*10^6;
fs_lw = sps*fbaud;
t_lw= t_duration;

%Generating from laser
%Pin to laser is the avg output power
Ein= LASER(Pin,lw,int64(fs_lw*t_lw), fs_lw);

%Passing to MZM
Eout_tx= MZM(Ein,Vpi,V_sig,V_off,IL,ER,Pin); 

%Transmitting through fiber
L=80*10^3;
Efib_out = fiber_out(L,fs_lw,Eout_tx);
Pfib_out = (abs(Efib_out)).^2;

%Input to Avalanche Photodetector
%PD params
PD_bw = 30*10^9; 
Rd = 0.9; 
RL = 50; 
%Defining known constts.
q = 1.6*10^(-19);
k = 1.38064852 *10^(-23);
T = 300;
Fn = 1; %Assuming NF from amplifiers is 0dB 
M = 1000; %Gain of APD
FA = M^0.7; %Typical APD NF

%Finding output voltage signal from APD
N_sig = length(Efib_out);
Vout_rx = zeros(N_sig,1);
%Averaging output over iterations
iter = 5;
for i=1:iter
Vout_rx =  Vout_rx+APD(Efib_out, Rd, PD_bw, fs_lw, RL,M,FA,att); %Noisy volt signal
end
Vout_rx = Vout_rx/iter;

%Constellation plot
%Choosing middle samples in every symbol slot
V_rx_constt = Vout_rx(floor(sps/2):sps:end);
V_rx_constt = V_rx_constt/(max(abs(V_rx_constt)));
scatterplot(V_rx_constt)
title('Voltage Constellation output with an APD of gain M = 1000 with RL = 50 ohms')
grid on

%Finding thoeretical avg. SNR for the PD
SNR =  ((M*Rd*mean(Pfib_out)))^2 /( 2*q*(M^2)*FA*mean(Rd*Pfib_out)*PD_bw +4*k*T*Fn*PD_bw/(RL) );
fprintf( 'SNR(in dB)'+string(10*log10(SNR))+'\n')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AVALANCHE PHOTODETECTOR TRANSFER FUNCTION

%Ein- column vector
%bw is one sided PD bw.
%NOTE- We assume the Ein is due to transmitted signal + ambient sources

function Vout = APD( Ein, Rd, bw, fs,RL,M,FA,att)

%Finding power
Pin = (abs(Ein)).^2;

%Defining noiseless o/p current (including dark current)
Id_1 = M*Rd*Pin;
I_av = mean(Id_1);

%Defining known constts.
q = 1.6*10^(-19);
k = 1.38064852 *10^(-23);
T = 300;
Fn = 1; %Assuming NF from amplifiers is 0dB 

%shot and thermal noise is given by:
%shot: sigma_s^2 =  2q*M^2*FA*Rd*Pin_tot* bw (of PD)
%Thermal: sigma_t^2 = 4kT*Fn*bw/RL 
%We also have to traverse through the whole signal

N = length ( Ein);
Id_net = zeros(N,1);
noise_curr = zeros(N,1);


%Traversing through samples and adding noise accordingly with variance
%having the delta_f of fs/N. N samples would give variance with delta f
%being fs. Now, we add a Low pass filter in the end to restrict it to the
%bandwidth of the PD (bw). Effectively making the delta_f overall to be bw of PD.
for i = 1 :N
%Net noise given factors M and Noise factor FA
%Note- the factor is M*FA to the amplified current Id =  M*Rd*Pin
sigma_t= sqrt( 2*q*M*FA*Id_1(i)*fs/N +4*k*T*Fn*fs/(N*RL));
%Mulitplying std dev required to a unit variance white noise
noise_curr(i) = sigma_t*randn;
%Adding white noise of given variance to signal 
Id_net(i) = Id_1(i) + noise_curr(i) ;
end


%Low pass filter with BW = PD one sided BW
Id_net_lpf = lowpass(Id_net,bw,fs,'StopbandAttenuation',att);

%Net voltage
Vout = Id_net_lpf*RL;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FIBER TRANSFER FUNCTION
function E_out_fib= fiber_out(L,fs_lw,Ein_fib)

%params of transfer function of dispersion
lamda = 1550*10^-9;
c = 3*10^8;

%Dispersion coeff typical value of SMF
D = 17*10^(-6); % converted from 17 ps/(nm-km) to SI


%calculating FFT of input sig.
%fs_lw is sampling freq of the E-field laser output

fft_inp = fft(Ein_fib);
fft_inp=fftshift(fft_inp);
%constructing transfer function
freq = linspace(-fs_lw/2, fs_lw/2,length(fft_inp)+1);
freq = freq(1:end-1);

%NLSE
%here omega is 2pi*f
omega = 2*pi*freq;
T = zeros(length(omega),1);
theta = zeros(length(omega),1);
for i=1:length(omega)
theta(i) = -D*(lamda^2)*((omega(i))^2)*L/(4*pi*c);
%T(i) = cos(theta)+1j*sin(theta);
T(i) = exp(1j*theta(i));
end

%output sig.
fft_Eout_fib = T.*fft_inp; 


%output time domain waveform
%ifft found after reverting the earlier fftshift in fourier domain
E_out_fib = ifft(ifftshift(fft_Eout_fib));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MZM FUNCTION
function E_out = MZM(E_in,Vpi,V,V_off,IL,ER,Pin)
%Pin is in W
%10*log10(Pout_1) = 10*log10(Pin)-IL;
Pout_1 = 10^(-IL/10) * Pin;
Pout_0 =  Pout_1* 10^(-ER/10);
%Enet = (Ein/2)*(1+ eta*e^(j theta))
% ER = (1+eta)^2 / ( 1-eta)^2
% ER^0.5 = (1 +eta)/ (1-eta) 
%(ER^0.5 -1)/(ER^0.5 +1) = eta 

ER1 = 10^(ER/10); %ratio form
eta = (ER1^0.5 -1)/(ER1^0.5 +1);
%theta/2 = V*pi/(2*Vpi)
%theta = V*pi/(Vpi)
theta = (V-V_off)*pi/Vpi;
%The following is done to get right BPSK mod.
Enet = (E_in/2).*(cos(-theta/2)+j*sin(-theta/2)) + eta*(E_in/2).*(cos(theta/2)+j*sin(theta/2));

%Insertion loss factor
%Enet max  =  (Ein/2) * (1 +eta)
Pnet_max = Pin*((1+eta)^2)/4;
Pout_1 = 10^(-IL/10) * Pin;
%in terms of power
f_iL = Pout_1/Pnet_max;
% factor for E field
f_iL = f_iL^0.5;

E_out = f_iL*Enet;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LASER function
function field = LASER(PAVG,LW,LEN,FS)
% Phase noise genc'
rand_var = randn(LEN,1);
sigma = sqrt(2*pi*LW.*(1/FS)); 
noise_vec = (ones(LEN,1) .* sigma) .* rand_var;      
noise_vec(1)=0;
phase_noise=cumsum(noise_vec,1); %
field = ((PAVG)^0.5).*(exp(1i*phase_noise)) ; 
end


