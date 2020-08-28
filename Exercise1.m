clc
clear all

N = [10,20,40,80,160,320]; 
T = cell(1,length(N));
U = cell(1,length(N));
f=@(t,u) t-(t*(u^2))-u;


for i = 1 : length(N)

    h = 1/N(i); %Step Length 
    [t, u] = RK(f, 1, h);
    T{i} = t;
    U{i} = u;
end

plot(T{5},U{5} );
