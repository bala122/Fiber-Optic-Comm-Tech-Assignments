
%Q1 
lw= 10^3;

%a)
Vpi = 3;
V_off = 0;
IL = 2;
ER = 30;
Pin = -10; %in dBm
%Converting to W
Pin = 10^(Pin/10);
Pin = Pin *10^(-3);

%V signal to MZM
% V_off + V_bias + V_modulation
%Here, we need only one bias point each iteration so,
% it'll be V_off + Vbias_i

%Vbiasi shd be -2*Vpi to 2*Vpi to include 3 peaks 2 nulls

Vb = linspace(-2*Vpi, 2*Vpi,400); % 400 points for enough resolution
V_sig = V_off + Vb;
P_out = zeros(1, length(Vb));
for i = 1 :length(Vb)
Eini= LASER(Pin,lw,1, 50*10^3);
E_outi= MZM(Eini,Vpi,V_sig(i),0,IL,ER,Pin);
P_out(i) = (abs(E_outi))^2;
end

%Intensity plotted is only |E|^2 , ignoring other constants involving area
%, etc.

plot(V_sig, P_out)
title('Plot of Intensity vs Voltage input for the MZM ( a))')
xlabel('Voltage (V) ->')
ylabel(' Intensity (V/m)^2 ->')
grid on


%b)
Vpi = 4;
V_off = 0.8;
IL = 3;
ER = 25;
Pin = -10; %in dBm
%Converting to W
Pin = 10^(Pin/10);
Pin = Pin *10^(-3);

%Vbiasi shd be -2*Vpi to 2*Vpi to include 3 peaks 2 nulls

Vb = linspace(-2*Vpi, 2*Vpi,400); % 400 points for enough resolution
V_sig = V_off + Vb;
P_out = zeros(1, length(Vb));
for i = 1 :length(Vb)
Eini= LASER(Pin,lw,1, 50*10^3);
E_outi= MZM(Eini,Vpi,V_sig(i),V_off,IL,ER,Pin);
P_out(i) = (abs(E_outi))^2;
end

%Intensity plotted is only |E|^2 , ignoring other constants involving area
%, etc.

figure
plot(V_sig, P_out)
title('Plot of Intensity vs Voltage input for the MZM ( b))')
xlabel('Voltage (V) ->')
ylabel(' Intensity (V/m)^2 ->')
grid on






function E_out = MZM(E_in,Vpi,V,V_off,IL,ER,Pin)

%For the extinction ratio, we assume some 'eta' factor to signify loss in
%the 2nd path. We then use this eta to relate to ER. 
%Then, for insertion loss, we just scale the whole plot / intensity values
%by the insertion loss factor so that even the maxima gets scaled acc. to
%insertion loss and the ER is still maintained.

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