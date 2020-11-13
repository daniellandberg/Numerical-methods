steps = 5;
[A, h] = getA(steps);
dt = 1*h^2;
tvec = (0:dt:2);
Temperatur = zeros(length(tvec),(steps+2));
h1=1/7;
x = [0:h1:1]; 
t = [0:dt:2];
%Euler
for i = 1 : length(tvec)-1
    b = (funcB(tvec(i), (steps+2)))/h^2;
    Temperatur(i+1,:) = Temperatur(i,:)' + dt*(A*Temperatur(i,:)'+b);
end

%add x0
x0 = zeros(length(tvec),1);
for i = 1 : length(tvec)
    y = funcB(tvec(i), (steps+2));
    x0(i) = y(1);
end
Matrix = zeros(length(tvec),(steps+3));
Matrix(:,1)=x0;
for i = 2 : (steps+3)
    Matrix(:,i) = Temperatur(:,i-1);
end
mesh(x,t,Matrix)
xlabel('X');
ylabel('time')
title('mesh plot showing the heat change over time')
figure(2)
TimeZerosPointFive = Matrix(46,:);
VecToPlot = (0:1:7);
plot(VecToPlot,TimeZerosPointFive);
xlabel('X');
ylabel('temperature')
title('tempeture difference at time 0.5')


figure(3)
firstx = Matrix(:,1);
lastx = Matrix(:,end);
plot(tvec,firstx)
hold on
plot(tvec,lastx)
title('tempeture difference at X0 and Xn')
legend('X0', 'Xn')
xlabel('time');
ylabel('temperature')
%%
%2b)

%u0 = zeros((steps+2),1);
tic
step = 160;
[A, h] = getA(step-2);
u = zeros(step,1);
R=@(t,u) A*u+gb(t, step);
options=odeset('Jacobian',A);
[t,y] = ode23s(R,[0 2],u,options);
toc
mesh(y)
%%
%3
T00 = 20;
T52 = 200;
N = 49; %9
M = 21; %3
h = 0.1;
x = [0:h:5-2*h]; 
y = [0:h:2];

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

for i = 1:N
    A(i, N+i) = -2; %changing A to match the neuman BC
end
dia = (N*M)+1-(2*N);
for i = ((N*M)-N)+1: N*M %changing A to match the neuman BC
    A(i, dia) = -2;
    dia = dia + 1;
end
A = A/h^2;

% Time vector
dt = 0.1;
tvec = (0:dt:10);
% temperature vector
Temperatur = zeros(length(tvec), N*M);

% IC
b0 = getu0(N,M, T00, T52, h);
u0 = A\b0; %from 2a in previous exercise
Temperatur(1,:) = u0;

% Crank-Nicolson
I = eye(N*M);
b = getB(N,M, T00, T52, h);
tempers = zeros(N*M,1); % u_t = 0
A = -A;
VAL=reshape(u0,N,M);
threeOne=zeros(length(tvec),1);
threeOne(1) = VAL(30,10);
for i = 1 : length(tvec)-1
    known = (I + ((dt*A)/2))*Temperatur(i,:)' + (dt*b); %from lecture 10
    Spars = sparse((I - ((dt/2)*A))); %from lecture 10
    Temperatur(i+1,:) = Spars\known; %from lecture 10
    
    Tempi=reshape(Temperatur(i+1,:),N,M);
    threeOne(i+1)=Tempi(30,10); 
    % Plot in every time-step
%     Tempi=reshape(Temperatur(i,:),N,M);
%     mesh(y,x,Tempi)
%     pause(0.1)
end
% 1 using delta x = delta y = delta t = 0.1
% Plotting t = 0, 0.5, 2, 10
figure(1)
Tempi=reshape(Temperatur(end,:),N,M);
mesh(y,x,Tempi)
title('t = 10')
xlabel('y');
ylabel('x')

% Plotting u(3,1,t) as a function of time
figure(2)
plot(tvec, threeOne);
title('u(3,1,t)')
xlabel('time')
ylabel('heat')


%%

function b = getB(N,M, T00, T52,h)
b = zeros((N*M),1);
F=@(x,y) 50+500*exp((-2*((x-1)^2))-((y-1.5)^2)); %eq

counter = 1;
for i = 1 : M %generate f using function F(x,y)
    for j = 1 : N
        b(counter) = F(j*h,i*h-h); 
        counter = counter +1;
    end
end


% Set boundary conditions when t > 0:
V = 1:N:length(b);
V2 = N:N:length(b);
b(V,1) = b(V,1)+(T00/(h^2));
b(V2,1) = b(V2,1)+(T52/(h^2));
%prepare to change f for the values where the Neuman BC applies
aj0 = zeros(N,1);
aj2 = zeros(N,1);
for i = 1 : N
    aj0(i) = F(i*h,0);
    aj2(i) = F(i*h,2);
end
Down = 1:1:N;
b(Down,1) = b(Down,1)-((-1)*2*h*aj0); % BC neuman y = 0
Up = (N*M-N)+1:1:(N*M);
b(Up,1) = b(Up,1)+((-1)*2*h*aj2); %BC neuman y = 2
end
function u0 = getu0(N,M, T00, T52,h)
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
u0 = f/h^2;


end

function [A, h] = getA(n)

L=1;
h = L/(n+1);

atemp = 1; %values for diag(a, -1)
aj=atemp*ones(1,n+1); 
btemp = -2; %values for diag(b)
bj=btemp*ones(1,n+2);
ctemp = 1; %values for diag(c, 1)
cj=ctemp*ones(1,n+1);

A1 = diag(bj);
A2 = diag(aj, -1);
A3 = diag(cj, 1);
A = (A1+A2+A3)/h^2; %the matrix with a, b and c on the diag
Areplace = zeros(length(bj),1);
Areplace(length(Areplace)-1) = 2/h^2;  %/h^2????
Areplace(end) = -2/h^2;  %/h^2????
A(end,:) = Areplace';
end

function b = funcB(vec, steps)
q=@(t) sin(pi*t/2);
b = zeros(steps,1);
if vec <= 1
    b(1) = q(vec);

end
end

function b = gb(t, step)
q=@(t) sin(pi*t/2);
b = zeros(step,1);
if t <= 1
    b(1) = q(t);
end
end
