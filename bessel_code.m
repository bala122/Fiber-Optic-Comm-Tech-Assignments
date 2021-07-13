b = linspace(0,1,10000);
V=7.5;

%l=0
plot(b,(V*((1-b).^0.5)).*besselj(1,V*((1-b).^0.5))./besselj(0,V*((1-b).^0.5)))
grid on
hold on
title('l=0')
plot(b,(V*((b).^0.5)).*besselk(1,V*((b).^0.5))./besselk(0,V*((b).^0.5)))
% ylim([-1000,1000])


figure

%l>=1
l=1
plot(b,(V*((1-b).^0.5)).*besselj(l-1,V*((1-b).^0.5))./besselj(l,V*((1-b).^0.5)))
grid on
hold on
title('l>=1')
plot(b,-(V*((b).^0.5)).*besselk(l-1,V*((b).^0.5))./besselk(l,V*((b).^0.5)))
% ylim([-1000,1000])