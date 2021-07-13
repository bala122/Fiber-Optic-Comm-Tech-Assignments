

%% PD parameters
PD_bw = 12*10^9; %Typical BW of PD = 12GHz (above fbaud)
Rd = 0.5; %Typical value of 0.5 A/W
RL = 50; %Typical load resistance
%Defining known constts.
q = 1.6*10^(-19);
k = 1.38064852 *10^(-23);
T = 300;
Fn = 1; %Assuming NF from amplifiers is 0dB 
att = 500; %attenuation of close to ideal LPF

%% GENERATION OF E-FIELD

%Generating PRBS

%QPSK parameters
%Input Data Parameters
%Data rate 20Gbps
%fbaud 10Gbaud
%sps = 4
%No. of symbols 2^12

no_symb = 2^12;
y_prbs9= prbs(9,no_symb);
y_prbs9_2 = circshift(y_prbs9,2);
f_baud = 10*10^9;
t_baud = 1/f_baud;
Nt1 = 4; %sps
fs = Nt1*f_baud;
% Repeating and then concatenating
y_p9 = repmat(transpose(y_prbs9),1,Nt1);
y_p9_2 = repmat(transpose(y_prbs9_2),1,Nt1);

prbs_seq=[];
prbs_seq_2=[];
for i=1:no_symb
    prbs_seq = [ prbs_seq y_p9(i,:)];
    prbs_seq_2 = [ prbs_seq_2 y_p9_2(i,:)];
end

%Power input

%The  " power " = abs(E)^2 calculated from the following value would be
%about 0dBm (peak) which sounds reasonable for a test power given typical
%NEP values of photodetectors .
P_IQ = 1*10^(-3);
E_field = ((P_IQ/2)^0.5 )*( 2.*(prbs_seq-0.5)+1j.*2*(prbs_seq_2-0.5) ) ;

%% TESTING AND DEMODULATION

%ELO
%Generating ELO with 10 times power
Elo = ((10*P_IQ)^0.5)*ones(1,Nt1*no_symb);
[i_I,i_Q] = coh_det(E_field,Elo,Rd,PD_bw,fs,att);

i_I_constt = i_I(floor(Nt1/2):Nt1:end);
i_Q_constt = i_Q(floor(Nt1/2):Nt1:end);

i_IQ_max = max(  max(abs(i_I_constt)) , max(abs(i_Q_constt)) ); 
scatterplot( (i_I_constt+1j.*i_Q_constt)/i_IQ_max )
title('Constellation diagram output of coherent detector')
grid

%% COHERENT DETECTOR
function [i_I,i_Q] = coh_det(Es,Elo,Rd,pd_bw,fs,att)
% Inputs- Es, ELO
% Outputs- E1, E2 (in-phase) E3, E4 (Q-phase)

%In phase 
E1 = 0.5*(Es+1j.*Elo);
E2 = 0.5*(1j.*Es+Elo);
i1 = PD(E1 ,Rd, pd_bw,fs,50,att); 
i2 = PD(E2 ,Rd, pd_bw,fs,50,att);
i_I = i1-i2;

%Q-phase
E3 = 0.5*(1j.*Es-1j.*Elo);
E4 = 0.5*(-Es-Elo);
i3 = PD(E3 ,Rd, pd_bw,fs,50,att);
i4 = PD(E4 ,Rd, pd_bw,fs,50,att);
i_Q = i4-i3;
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

%Low pass filter with BW = PD one sided BW
Id_net_lpf = lowpass(Id_net,bw,fs,'StopbandAttenuation',att);

%Net current
Iout = Id_net_lpf;

end
