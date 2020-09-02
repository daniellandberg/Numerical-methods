clc
clear all

N = [10,20,40,80,160,320]; 
T = cell(1,length(N));
U = cell(1,length(N));
f=@(t, y) [y(2); (1-y(1)^2)*y(2)-y(1)];


for i = 1 : length(N)

    h = 1/N(i); %Step Length 
    [t, u] = RK(f, 1, h);
    T{i} = t;
    U{i} = u;
end
A = U{1}(1,:);
t = [t 1];
plot(t, A);
