function [t, y] = RK(f, y0, inVals,N)

 %vector with 320 even steps between 0 and 1
y=zeros(inVals,N+1);
y(:,1)=y0;
h=1/N;
t=linspace(0,1,N+1)';




for i = 1:N
    k1 = f(y(:,i));
    k2 = f(y(:,i)+h*k1);
    k3= f(y(:,i)+(h*k1/4)+(h*k2/4));
    y(:,i+1) = y(:,i) + ((h/6)*(k1+k2+k3));
    
end

