clc
clear all
Z = cell(1,4);
U = cell(1,4);
n = [9, 19, 39, 79];
%n = [39,79, 159,319];
v = [0.1, 0.5,1 , 10];
%n = [159,319, 639,1279,2559]; %n for order of accuracy
for i = 1 : length(n)
    [z, u]=fda(n(i), 1);
    Z{i} = z;
    U{i} = u;   
end
figure(1)
plot(Z{1}, U{1})
hold on
plot(Z{2}, U{2})
hold on
plot(Z{3}, U{3})
hold on
plot(Z{4}, U{4})

legend('h = 1', 'h = 0.5', 'h = 0.25', 'h = 0.125')
title('V=1')
xlabel('z');
ylabel('T')

%%
en = zeros(1,length(n)-1);
tend = u(end);
for i = 1 : length(n)-1
    temp = U{i}(:,1);
    yn = temp(end);
    en(i) = abs(tend - yn);
end

hLista = ones(1,length(n)-1);
for i = 1 : length(n)-1
   hLista(i) = hLista(i)/(n(i)+1); 
end
figure(2);
loglog(hLista,en);
%%

T00 = 20;
T52 = 200;
N = 49;
M = 21;
h = 0.1;
A = zeros(N*M);
for i = 1:N*M
    A(i,i) = 4;
end

for i = 1:N-1
    for j = 1:M
        A(i+(j-1)*N, i+(j-1)*N+1) = -1; %around the 4
        A(i+(j-1)*N+1, i+(j-1)*N) = -1;
    end
end

for i = 1:N
    for j = 1:M-1
        A(i+(j-1)*N, i+j*N) = -1; %N steps away
        A(i+j*N, i+(j-1)*N) = -1;
    end
end

f = zeros((N*M),1);
% Set boundary conditions:
V = 1:N:length(f);
V2 = N:N:length(f);
f(V,1) = f(V,1)+T00;
f(V2,1) = f(V2,1)+T52;
Down = 1:1:N;
f(Down,1) = f(Down,1)-((-1)*2*h*0); %-2haj times the row where y = 0
Up = (N*M-N):1:(N*M);
f(Up,1) = f(Up,1)+((-1)*2*h*0); %2haj times the row where y = 2 
for i = 1:N
    A(i, N+i) = -2;
end
dia = (N*M)+1-(2*N);
for i = ((N*M)-N)+1: N*M
    A(i, dia) = -2;
    dia = dia +1;
end
A = A*1/h^2;
f = f/(h^2);
u = A\f;
u=reshape(u, N,M);
figure(3)
imagesc(u')
%%

h = [0.2,0.1, 0.05,0.025];
n = [24,49,99,199];
m = [11,21,41, 81,161];
U = cell(1,4);
for i = 1 : length(h)
    u = nog(n(i), m(i), h(i));
    U{i}=u;
end

%%
en = zeros(1,length(n)-1);
tend = u(end);
for i = 1 : length(n)-1
    temp = U{i}(:,1);
    yn = temp(end);
    en(i) = abs(tend - yn);
end

hLista = ones(1,length(n)-1);
for i = 1 : length(n)-1
   hLista(i) = hLista(i)/(n(i)+1); 
end
figure(4);
loglog(hLista,en);
%%

T00 = 20;
T52 = 200;
N = 49;
M = 21;
h = 0.1;
A = zeros(N*M);
for i = 1:N*M
    A(i,i) = 4;
end

for i = 1:N-1
    for j = 1:M
        A(i+(j-1)*N, i+(j-1)*N+1) = -1; %around the 4
        A(i+(j-1)*N+1, i+(j-1)*N) = -1;
    end
end

for i = 1:N
    for j = 1:M-1
        A(i+(j-1)*N, i+j*N) = -1; %N steps away
        A(i+j*N, i+(j-1)*N) = -1;
    end
end

f = zeros((N*M),1);
temp=zeros(N,M);
F=@(x,y) 50+500*exp((-2*((x-1)^2))-((y-1.5)^2)); %eq

counter = 1;
for i = 1 : M %generate f using function F(x,y)
    for j = 1 : N
        f(counter) = F(j*h,i*h-h); 
        counter = counter +1;
    end
end

% Set boundary conditions:
V = 1:N:length(f);
V2 = N:N:length(f);
f(V,1) = f(V,1)+(T00/(h^2));
f(V2,1) = f(V2,1)+(T52/(h^2));
%prepare to change f for the values where the Neuman BC applies
aj0 = zeros(N,1);
aj2 = zeros(N,1);
for i = 1 : N
    aj0(i) = F(i*h,0);
    aj2(i) = F(i*h,2);
end
Down = 1:1:N;
f(Down,1) = f(Down,1)-((-1)*2*h*aj0); % BC neuman y = 0
Up = (N*M-N)+1:1:(N*M);
f(Up,1) = f(Up,1)+((-1)*2*h*aj2); %BC neuman y = 2

for i = 1:N
    A(i, N+i) = -2; %changing A to match the neuman BC
end
dia = (N*M)+1-(2*N);
for i = ((N*M)-N)+1: N*M %changing A to match the neuman BC
    A(i, dia) = -2;
    dia = dia +1;
end
A = A*1/h^2;
u = A\f;
u=reshape(u, N,M);
figure(5)
imagesc(u')
%%
function [z,u] = fda(n, V)
v=V;
L = 10;
a = 2;
b = 4;
Q0 = 10;
k = 0.5;
alpha0 = 8;
ro = 1;
C = 0.75;
Tout = 20;
T0 = 50;

alphav = ((((v^2*ro^2*C^2)/4)+alpha0^2)^0.5)-((v*ro*C)/2);
h = L/(n+1);

atemp = -((k/h^2)+((v*ro*C/(2*h)))); %values for diag(a, -1)
aj=atemp*ones(1,n+1); 
btemp = (2*k)/h^2; %values for diag(b)
bj=btemp*ones(1,n+2);
ctemp = (-k/h^2)+(v*ro*C/(2*h)); %values for diag(c, 1)
cj=ctemp*ones(1,n+1);

A1 = diag(bj);
A2 = diag(aj, -1);
A3 = diag(cj, 1);
A = A1+A2+A3; %the matrix with a, b and c on the diag
Areplace = zeros(length(bj),1);
Areplace(length(Areplace)-1) = atemp+ctemp*1;
Areplace(end) = btemp - (ctemp*(2*h*alphav/k));
A(end,:) = Areplace';
A = sparse(A);

z=linspace(0,L,(n+2))';
q=@(z) Q0*sin(((z-a)*pi)/(b-a));
indexa = (a/h)+1; %find the index for the coil in step vector z
indexb = (b/h)+1;
Q=zeros(length(z),1);

Q(indexa:indexb) = q(z(indexa:indexb)); % put in the values of q(z) on index a/b: b/h in Q


Q(1) = Q(1) - (atemp*T0);
d1 = 2*h*alphav*(Tout/k);
Q(end) = Q(end) - (ctemp*d1); %calcuate boundary for Tout using the point Cn* alphav from BC equation

u = A\Q; %Solving the system
u(1) = T0;
%u = [T0;u]; % adding initial values?
z=linspace(0,L,(n+2))';
end

function [u] = nog(n, m, h)
T00 = 20;
T52 = 200;
N = n;
M = m;
h
A = zeros(N*M);
for i = 1:N*M
    A(i,i) = 4;
end

for i = 1:N-1
    for j = 1:M
        A(i+(j-1)*N, i+(j-1)*N+1) = -1; %around the 4
        A(i+(j-1)*N+1, i+(j-1)*N) = -1;
    end
end

for i = 1:N
    for j = 1:M-1
        A(i+(j-1)*N, i+j*N) = -1; %N steps away
        A(i+j*N, i+(j-1)*N) = -1;
    end
end

f = zeros((N*M),1);
% Set boundary conditions:
V = 1:N:length(f);
V2 = N:N:length(f);
f(V,1) = f(V,1)+T00;
f(V2,1) = f(V2,1)+T52;
Down = 1:1:N;
f(Down,1) = f(Down,1)-((-1)*2*h*0); %-2haj times the row where y = 0
Up = (N*M-N):1:(N*M);
f(Up,1) = f(Up,1)+((-1)*2*h*0); %2haj times the row where y = 2 
for i = 1:N
    A(i, N+i) = -2;
end
dia = (N*M)+1-(2*N);
for i = ((N*M)-N)+1: N*M
    A(i, dia) = -2;
    dia = dia +1;
end
A = A*1/h^2;
A = sparse(A);
f = f/(h^2);
u = A\f;
u=reshape(u, N,M);
end