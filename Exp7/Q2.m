

%% PD parameters
PD_bw = 26*10^9; %Typical BW of PD = 26GHz (above fbaud)
Rd = 0.5; %Typical value of 0.5 A/W
RL = 50; %Typical load resistance
%Defining known constts.
q = 1.6*10^(-19);
k = 1.38064852 *10^(-23);
T = 300;
Fn = 1; %Assuming NF from amplifiers is 0dB 
att = 500; %attenuation of close to ideal LPF

%% GENERATION OF QPSK E-FIELD


%Generating PRBS

%QPSK parameters
%Input Data Parameters
%Data rate 50Gbps
%fbaud 25Gbaud
%sps = 10
%No. of symbols 2^12

no_symb = 2^12;
y_prbs9= prbs(9,no_symb);
y_prbs9_2 = circshift(y_prbs9,2);
f_baud = 25*10^9;
t_baud = 1/f_baud;
Nt1 = 10; %sps
fs = Nt1*f_baud;
% Repeating and then concatenating
y_p9 = repmat(transpose(y_prbs9),1,Nt1);
y_p9_2 = repmat(transpose(y_prbs9_2),1,Nt1);

V_sig=[];
V_sig2=[];
for i=1:no_symb
    V_sig = [ V_sig y_p9(i,:)];
    V_sig2 = [ V_sig2 y_p9_2(i,:)];
end

%MZM data
Vpi = 3;
V_off = 0.5;
IL = 0.3;
ER = 30;

%Choosing a swing of 2Vpi (peak-peak) - individual BPSK mod.
V_sig = 2*Vpi*(V_sig - 0.5);
V_sig2 = 2*Vpi*(V_sig2 - 0.5);

%Biasing @ Vpi
V_sig = V_off+Vpi + V_sig;
V_sig2 = V_off+Vpi + V_sig2;


V_sig = transpose(V_sig); % Ein, V_sig are column vectors
V_sig2 = transpose(V_sig2); % Ein, V_sig2 are column vectors

%Laser linewidth-0Hz for no phase noise
lw = 0;

Pin = 0.1*10^(-3);
Ein_MZM = LASER(Pin,lw,Nt1*no_symb,fs);
E_outI = MZM(Ein_MZM,Vpi,V_sig,V_off,IL,ER,Pin);
E_outQ = MZM(Ein_MZM,Vpi,V_sig2,V_off,IL,ER,Pin);
%Factor of 1/(2^0.5) added for power conservation
E_out_IQ = (1/(2^0.5))*(E_outI +1j.*E_outQ);

%Adding amplitude noise to whole spectrum according to OSNR (measured in 12.5GHz band)
OSNR = 30;
P_noise = (fs/(12.5*10^9))*(mean(abs(E_out_IQ))^2)/(10^(OSNR/10));
amp_noise = sqrt(P_noise)*randn(Nt1*no_symb,1);
E_IQ_noisy = (abs(E_out_IQ) + amp_noise).*(exp(1j*(angle(E_out_IQ))));

%% TESTING AND DEMODULATION

%ELO
PLo = 10^(-3);
Elo = LASER(PLo,lw,Nt1*no_symb,fs);
[i_I,i_Q] = coh_det(E_IQ_noisy,Elo,Rd,PD_bw,fs,att);
i_IQ = i_I+1j*i_Q;
i_I_constt = i_I(floor(Nt1/2):Nt1:end);
i_Q_constt = i_Q(floor(Nt1/2):Nt1:end);

i_IQ_max = max(  max(abs(i_I_constt)) , max(abs(i_Q_constt)) ); 
scatterplot( (i_I_constt+1j.*i_Q_constt)/i_IQ_max )
title('Constellation diagram for I,Q currents -QPSK 25Gbaud,OSNR 30dB ')
grid

%Deciding symbols via thresholding
IQ_demod = zeros(no_symb,1);
for i=1:no_symb
if(i_I_constt(i)>0)
    if(i_Q_constt(i)>0)
        IQ_demod(i) = 1+1j;
    else
        IQ_demod(i) = 1-1j;
    end
else
    if(i_Q_constt(i)>0)
        IQ_demod(i) = -1+1j;
    else
        IQ_demod(i) = -1-1j;
    end
end
end

%Finding BER
I_inp = V_sig - Vpi- V_off;
Q_inp = V_sig2 - Vpi -V_off;

%Accounting for inversion in MZM
I_inp = -(I_inp)/Vpi;
Q_inp = -(Q_inp)/Vpi;
I_inp = I_inp(floor(Nt1/2):Nt1:end);
Q_inp = Q_inp(floor(Nt1/2):Nt1:end);
IQ_inp = I_inp+1j*Q_inp;

err_bits=0;
for i =1:no_symb
if(real(IQ_inp(i)-IQ_demod(i))~=0)
    err_bits= err_bits+1;
end

if(imag(IQ_inp(i)-IQ_demod(i))~=0)
    err_bits= err_bits+1;
end

end
BER = err_bits/no_symb;
fprintf('BER(OSNR 30dB):'+string(BER)'+'\n')

%Comparison plots TX, demodulated symbols
figure
subplot(2,1,1)
stem(real(IQ_inp(1:10)))
ylabel('Input symbol')
title('In-phase part of Symbol (Tx and demodulated)(OSNR 30dB)')
hold 
subplot(2,1,2)
stem(real(IQ_demod(1:10)))
ylabel('Output demodulated symbol')
xlabel('Symbol index')

figure
subplot(2,1,1)
stem(imag(IQ_inp(1:10)))
ylabel('Input symbol')
title('Q-phase part of Symbol (Tx and demodulated)(OSNR 30dB)')
hold 
subplot(2,1,2)
stem(imag(IQ_demod(1:10)))
ylabel('Output demodulated symbol')
xlabel('Symbol index')

%% COHERENT DETECTOR
function [i_I,i_Q] = coh_det(Es,Elo,Rd,pd_bw,fs,att)
% Inputs- Es, ELO
% Outputs- E1, E2 (in-phase) E3, E4 (Q-phase)

%In phase 
E1 = 0.5*(Es+1j.*Elo);
E2 = 0.5*(1j.*Es+Elo);
i1 = PD(E1 ,Rd, pd_bw,fs,50,att); 
i2 = PD(E2 ,Rd, pd_bw,fs,50,att);
i_Q = i1-i2;

%Q-phase
E3 = 0.5*(1j.*Es-1j.*Elo);
E4 = 0.5*(-Es-Elo);
i3 = PD(E3 ,Rd, pd_bw,fs,50,att);
i4 = PD(E4 ,Rd, pd_bw,fs,50,att);
i_I = i4-i3;

end

%% PHOTODETECTOR 
function Iout = PD( Ein, Rd, bw, fs,RL,att)

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
%We also have to traverse through the whole signal

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

Id_net_phase = angle(Id_net);
%Low pass filter with BW = PD one sided BW
Id_net_lpf = lowpass(Id_net,bw,fs,'StopbandAttenuation',att);
Id_net_lpf = abs(Id_net_lpf).*(exp(1j*(Id_net_phase)));
%Net current
Iout = Id_net_lpf;

end

%% MACH ZENDER MODULATOR

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
%The following is done to get right QPSK, BPSK mod.
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

%% LASER
function field = LASER(PAVG,LW,LEN,FS)
% Phase noise genc'
rand_var = randn(LEN,1);
sigma = sqrt(2*pi*LW.*(1/FS)); 
noise_vec = (ones(LEN,1) .* sigma) .* rand_var;      
noise_vec(1)=0;
phase_noise=cumsum(noise_vec,1); %
field = ((PAVG)^0.5).*(exp(1i*phase_noise))   ; 
end