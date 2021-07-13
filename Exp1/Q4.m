%Q4
%Eg calculated is (in J):
Eg = 1.053703703703704 * 1.6 *10^-19 ;
%temperature
T = 300;
%boltzmann constt.
k = 1.38064852 *10^-23;
%plancks constt.
h = 6.62607004 * 10^-34;

e = exp(1);
%R = R0 (hf- Eg)^0.5 * exp( -( hf - Eg)/(kT))
%Choosing around 10^4 Thz as the resolution for sufficient accuracy
%Plotting only from lowest possible frequency- Eg/h. Frequencies below this
%would have 0 Rsp anyway.
f = linspace( Eg/h, (Eg +7*k*T)/h, 10^5); 
R = zeros(length(f),1);
for i= 1:length(f)
R(i) = ((h*f(i) - Eg)^0.5)* exp(- (h*f(i) - Eg)/(k*T));
end
%Normalizing
R = R/(max(R));
%plotting in Thz freq.
plot(10^-12*f,R)
grid on
title('Plot of Rate of Spont.Emission vs. freq. (@ T = 300K)')
xlabel('frequency,- >( THz)')
ylabel('Normalized Rsp ->')
legend('Rsp @ 300K')
%Finding 3dB BW
f1=0;
f2=0;
%Checking for frequencies close to 3dB point (error in Rsp_norm. <= 1e-3 for
% sufficient accuracy)
for i = 1:length(f)
    
    if( (abs(R(i)-0.5)<=1e-3) && (f(i) < (Eg +k*T/2)/h))
       f1 =f(i);
    end
    if( (abs(R(i)-0.5)<=1e-3) && (f(i) > (Eg +k*T/2)/h))
        f2 = f(i);
    end
end
f1
f2
f2-f1
%This is almost 11.25 THz 


