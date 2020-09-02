function [t, u] = RK(f,y0,h)

 %vector with 320 even steps between 0 and 1
t = linspace(0,1,100);
y(1) = 1;
y(2) = 0;
y=y';

for i = 1:100
    k1 = f(t(i), y(:,i));
    k2 = f(t(i)+h, y(:,i)+h*k1);
    
    k3= f((t(i))+(h/2), y(:,i)+(h*k1/4)+(h*k2/4));
    y(:,i+1) = y(:,i) + ((h/6)*(k1+k2+k3));
   
    
end
u = y;