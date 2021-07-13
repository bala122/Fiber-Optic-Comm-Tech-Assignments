%Q5
%NOTE- phase generation is a random process. The chosen values are such
%that the phase spread is within a quadrant for one case and becomes a full
%circle for another case for most of the runs, although not all. 
%So, running it multiple times would give a better picture of the phase
%spread.


%MZM data
Vpi = 3;
V_off = 0.5;
IL = 2;
ER = 25;
Pin = -10; %in dBm
%Converting to W
Pin = 10^(Pin/10);
Pin = Pin *10^(-3);

%%%%%%%%%%%%%%%%%%%%%
%Time duration 50ns.
%sps chosen for this time duration.

sps = 4;
%data sequence for the modulation signal for each BPSK part
%1Gbaud QPSK implies 1Gbps for each of the BPSK channels
% PRBS13 sps=4 Bitrate 1Gbps
%Hence bit slot time is 1ns
%Sampling frequency is fixed as sps*bitrate = 4GHz.


%LW 1MHz
%Generate the Time domain modulating waveform for 50 symbols =( 50ns)
t_duration1 = 50*10^(-9);
no_symb= 50;

y_prbs13= prbs(13,no_symb);
%circhifted to not produce same identical sequences.
y_prbs13_2 = circshift(y_prbs13,4);
y_p13 = repmat(transpose(y_prbs13),1,sps);
y_p13_2 = repmat(transpose(y_prbs13_2),1,sps);

%Repeating symbols to generate time domain waveform.
V_sig=[];
V_sig2=[];
for i=1:no_symb
    V_sig = [ V_sig y_p13(i,:)];
    V_sig2 = [ V_sig2 y_p13_2(i,:)];
end

%Choosing a swing of 2Vpi (peak-peak) - individual BPSK mod.
V_sig = 2*Vpi*(V_sig - 0.5);
V_sig2 = 2*Vpi*(V_sig2 - 0.5);

%Biasing @ Vpi
V_sig = V_off+Vpi + V_sig;
V_sig2 = V_off+Vpi + V_sig2;


V_sig = transpose(V_sig); % Ein, V_sig are column vectors
V_sig2 = transpose(V_sig2); % Ein, V_sig2 are column vectors


%Laser data 
lw1 = 10^6;
lw2 = 10*10^6;
fs_lw = sps*10^9; %both lw's have same fs.
t_lw1= t_duration1;
t_lw2= t_duration1;

%Using one input electric field and passing it to two MZMs which have two
%different input voltage signal waveforms.

Ein_lw1= LASER(Pin,lw1,int64(fs_lw*t_lw1), fs_lw);
E_out1= MZM(Ein_lw1,Vpi,V_sig,V_off,IL,ER,Pin);
E_out2= MZM(Ein_lw1,Vpi,V_sig2,V_off,IL,ER,Pin);
%Adding a phase of 90 deg to output of one of the MZMs.
E_out_net_lw1 = E_out1 + j*E_out2;
%Picking points for constellation diag.
E_constt = E_out_net_lw1(int64(sps/2):sps:end);
scatterplot(E_constt/max(abs(E_constt)))
title('Constellation diagram QPSK LW 1MHz time- 50ns')
grid on


Ein_lw2= LASER(Pin,lw2,int64(fs_lw*t_lw2), fs_lw);
E_out1= MZM(Ein_lw2,Vpi,V_sig,V_off,IL,ER,Pin);
E_out2= MZM(Ein_lw2,Vpi,V_sig2,V_off,IL,ER,Pin);
E_out_net_lw2 = E_out1 + j*E_out2;
%Picking points for constellation diag.
E_constt = E_out_net_lw2(int64(sps/2):sps:end);
scatterplot(E_constt/max(abs(E_constt)))
title('Constellation diagram QPSK LW 10MHz time- 50ns')
grid on


%TIME DURATION 600ns %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate the Time domain modulating waveform for 600 symbols =( 600ns)
t_duration2 = 600*10^(-9);
no_symb= 600; %no. of symbols
sps = 8;

y_prbs13= prbs(13,no_symb);
%circhifted to not produce same identical sequences.
y_prbs13_2 = circshift(y_prbs13,4);
y_p13 = repmat(transpose(y_prbs13),1,sps);
y_p13_2 = repmat(transpose(y_prbs13_2),1,sps);

V_sig=[];
V_sig2=[];
for i=1:no_symb
    V_sig = [ V_sig y_p13(i,:)];
    V_sig2 = [ V_sig2 y_p13_2(i,:)];
end

%Choosing a swing of 2Vpi (peak-peak)
V_sig = 2*Vpi*(V_sig - 0.5);
V_sig2 = 2*Vpi*(V_sig2 - 0.5);

%Biasing @ Vpi
V_sig = V_off+Vpi + V_sig;
V_sig2 = V_off+Vpi + V_sig2;


V_sig = transpose(V_sig); % Ein, V_sig are column vectors
V_sig2 = transpose(V_sig2); % Ein, V_sig2 are column vectors


%Laser data 
lw1 = 10^6;
lw2 = 10*10^6;
fs_lw = sps*10^9; 
%fs determined by sps is different for different constellation spreads
%for different times
t_lw1= t_duration2;
t_lw2= t_duration2;



Ein_lw1= LASER(Pin,lw1,int64(fs_lw*t_lw1), fs_lw);
E_out1= MZM(Ein_lw1,Vpi,V_sig,V_off,IL,ER,Pin);
E_out2= MZM(Ein_lw1,Vpi,V_sig2,V_off,IL,ER,Pin);
E_out_net_lw1 = E_out1 + j*E_out2;

%Picking constellation symbol points
E_constt = E_out_net_lw1(4:8:end);
scatterplot(E_constt/max(abs(E_constt)))
title('Constellation diagram QPSK LW 1MHz time- 600ns')
grid on


Ein_lw2= LASER(Pin,lw2,int64(fs_lw*t_lw2), fs_lw);
E_out1= MZM(Ein_lw2,Vpi,V_sig,V_off,IL,ER,Pin);
E_out2= MZM(Ein_lw2,Vpi,V_sig2,V_off,IL,ER,Pin);
E_out_net_lw2 = E_out1 + j*E_out2;

%Picking constellation symbol points
E_constt = E_out_net_lw2(2:4:end); 
scatterplot(E_constt/max(abs(E_constt)))
title('Constellation diagram QPSK LW 10MHz time- 600ns')
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


