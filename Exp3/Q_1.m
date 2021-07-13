%duration times
t1=10*10^-6;
t2=100*10^-6;
P_avg = 10^-3;
lw1= 10^3;
lw2= 10*10^3;
lw3= 100*10^3;
%We can adopt a fixed sampling rate for both times.
%min. points = 1000 implies fsamp*10 us = 1000 atleast.
%Hence sampling freq is:
fs = 100*10^6;

%Making E field output 

%LW 1KHz
%time 10us
E1_lw1= LASER(P_avg,lw1,int64(fs*t1), fs);
%Scatterplot on the E-field output
scatterplot(E1_lw1)
grid on
title('Plot of E-field run for 10us (LW = 1kHz)')
xlabel('Real part ->(V/m)')
ylabel('Imaginary part -> (V/m)')
legend('time: 10us LW = 1kHz')

%time 100us
E2_lw1= LASER(P_avg,lw1,int64(fs*t2), fs);
scatterplot(E2_lw1)
grid on
title('Plot of E-field run for 100us (LW = 1kHz)')
xlabel('Real part ->(V/m)')
ylabel('Imaginary part -> (V/m)')
legend('time: 100us LW = 1kHz')

%LW 10KHz
%time 10us
E1_lw2= LASER(P_avg,lw2,int64(fs*t1), fs);
scatterplot(E1_lw2)
grid on
title('Plot of E-field run for 10us (LW = 10kHz)')
xlabel('Real part ->(V/m)')
ylabel('Imaginary part -> (V/m)')
legend('time: 10us LW = 10kHz')

%time 100us
E2_lw2= LASER(P_avg,lw2,int64(fs*t2), fs);
scatterplot(E2_lw2)
grid on
title('Plot of E-field run for 100us (LW = 10kHz)')
xlabel('Real part ->(V/m)')
ylabel('Imaginary part -> (V/m)')
legend('time: 100us LW = 10kHz')


%LW 100KHz
%time 10us
E1_lw3= LASER(P_avg,lw3,int64(fs*t1), fs);
scatterplot(E1_lw3)
grid on
title('Plot of E-field run for 10us (LW = 100kHz)')
xlabel('Real part ->(V/m)')
ylabel('Imaginary part -> (V/m)')
legend('time: 10us LW = 100kHz')

%time 100us
E2_lw3= LASER(P_avg,lw3,int64(fs*t2), fs);
scatterplot(E2_lw3)
grid on
title('Plot of E-field run for 100us (LW = 100kHz)')
xlabel('Real part ->(V/m)')
ylabel('Imaginary part -> (V/m)')
legend('time: 100us LW = 100kHz')


function field = LASER(PAVG,LW,LEN,FS)
% Phase noise gen
rand_var = randn(LEN,1);
sigma = sqrt(2*pi*LW.*(1/FS)); 
noise_vec = (ones(LEN,1) .* sigma) .* rand_var;      
noise_vec(1)=0;
phase_noise=cumsum(noise_vec,1); %
field = ((PAVG)^0.5).*(exp(1i*phase_noise))   ; 
end