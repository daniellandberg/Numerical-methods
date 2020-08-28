function [t, u] = RK(f,y0,h)

t = 0:h:1; %vector with 320 even steps between 0 and 1
u(1) = y0;  %y(0)=1



for i = 1:length(t)-1
    k1 = f(t(i), u(i));
    k2 = f(t(i)+h, u(i)+h*k1);
    k3= f(t(i)+(h/2), u(i)+(h*k1/4)+(h*k2/4));
    u(i+1) = u(i) + ((h/6)*(k1+k2+k3));
    t(i+1) = (i+1)*h;
    
end