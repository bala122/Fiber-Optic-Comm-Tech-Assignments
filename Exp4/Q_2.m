%Q2

%MZM data
Vpi = 3;
V_off = 0.5;
IL = 2;
ER = 25;
Pin = -10; %in dBm
%Converting to W
Pin = 10^(Pin/10);
Pin = Pin *10^(-3);

%data sequence 
% PRBS9 sps=8 Bitrate 10Gbps
%Hence bit slot time is 0.1ns
%Sampling frequency is fixed as sps*bitrate = 80GHz.

%Generate the Time domain modulating waveform for 100 bits =( 10ns)
t_duration = 10*10^(-9);
no_symb= 100;


%Making the time domain input OOK voltage signal by
%repeating the symbol and making it look like a square wave.
y_prbs9= prbs(9,no_symb);
y_p9 = repmat(transpose(y_prbs9),1,8);
V_sig=[];
for i=1:no_symb
    V_sig = [ V_sig y_p9(i,:)];   
end

%Choosing a swing of Vpi (peak-peak)
V_sig = Vpi*(V_sig- 0.5);

%Biasing @ 0.5 Vpi
V_sig = V_off+Vpi/2+V_sig;

V_sig = transpose(V_sig); % Ein, V_sig are column vectors


%Laser data
lw= 100*10^3;
fs_lw = 80*10^9;
t_lw= t_duration;

Ein= LASER(Pin,lw,int64(fs_lw*t_lw), fs_lw);
E_out= MZM(Ein,Vpi,V_sig,V_off,IL,ER,Pin);

%Intensity plotted is only |E|^2 , ignoring other constants involving area
%, etc.

P_out = (abs(E_out)).^2;

phase_out = (180/pi).*(angle(E_out)); %in degrees
t1 = linspace( 0, t_duration, no_symb*8);

%subplot power and phase plots
subplot(2,1,1);
plot(10^9*t1,P_out)
title('Plot of output intensity vs time (OOK)')
xlabel('time (ns) ->')
ylabel('Intensity (normalized) (V/m)^2')
grid on

subplot(2,1,2);
plot(10^9*t1,phase_out)
title('Plot of output phase vs time (OOK)')
xlabel('time (ns) ->')
ylabel('Phase (deg.) ->')
grid on

%Constellation diagram
%Picking only the middle value of E-field in a bit slot for symbol
%constellation
E_constt = E_out(4:8:end);
%normalized.
scatterplot(E_constt/max(abs(E_constt)))
title('Constellation diagram OOK')
grid on



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
%LASER func
function field = LASER(PAVG,LW,LEN,FS)
% Phase noise genc'
rand_var = randn(LEN,1);
sigma = sqrt(2*pi*LW.*(1/FS)); 
noise_vec = (ones(LEN,1) .* sigma) .* rand_var;      
noise_vec(1)=0;
phase_noise=cumsum(noise_vec,1); %
field = ((PAVG)^0.5).*(exp(1i*phase_noise))   ; 
end




