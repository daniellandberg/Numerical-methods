
% assignment 1
N = [10,20,40,80,160,320]; 
T = cell(1,length(N));
U = cell(1,length(N));

f=@(y) [y(2); (1-y(1)^2)*y(2)-y(1)];

for i = 1 : length(N)
    [t, u] = RK(f, [1 0], 2, N(i), 1);
    T{i} = t;
    U{i} = u;   
end

%Order of accuracy
en = zeros(1,5);
temp = U{6}(1,:);
y1 = temp(end);
for i = 1 : length(N)-1
    temp = U{i}(1,:);
    yn = temp(end);
    en(i) = abs(y1 - yn);
end

hLista = ones(1,5);
for i = 1 : length(N)-1
   hLista(i) = hLista(i)/N(i); 
end
figure(1);
loglog(hLista,en);
title('Error plot of Runge kutta method')
xlabel('StepLength (h)')
ylabel('Error en')


%% 

%ex2.1
r1 = 0.04;
r2 = 10^4;
r3 = 3*10^7;
f = @(x)[-r1*x(1)+(r2*x(2)*x(3)); (r1*x(1))- (r2*x(2)*x(3))-(r3*(x(2)^2)); r3*(x(2)^2)];
N = [125, 250, 500, 1000, 2000];
Steg = [500, 510, 520, 530, 540, 550, 560, 570, 580, 580, 600]; %test
x0 = [1 0 0];

for i = 1 : length(N)
   [t, u] = RK(f, x0, 3,N(i), 1);
    T2{i} = t;
    U2{i} = u;
end

u  = U2{5};
t = T2{5};
setN = 5;

figure(4);
subplot(2,2,1);
loglog(T2{setN},U2{setN}(1,:)');
title('x1 when N = 1000')
xlabel('t');
ylabel('x1')


subplot(2,2,2);
loglog(T2{setN},U2{setN}(2,:)');
title('x2 when N = 1000')
xlabel('t');
ylabel('x2')

subplot(2,2,3);
loglog(T2{setN},U2{setN}(3,:)');
title('x3 when N = 1000')
xlabel('t');
ylabel('x3')

en = U2{5}(:,end)-U2{4}(:,end);

%%
%%2.b
x0 = [1 0 0];
r1 = 0.04;
r2 = 10^4;
r3 = 3*10^7;
options = odeset('RelTol', 10^-6);
f = @(t, x)[-r1*x(1)+(r2*x(2)*x(3)); (r1*x(1))- (r2*x(2)*x(3))-(r3*(x(2)^2)); r3*(x(2)^2)];
[t, x] = ode23s(f, [0 1], x0', options);


h=0;
for i = 1 : length(t)-1
    h(i) = abs(t(i) - t(i+1));
end

figure(11)
h = [h 0];
plot(t, h)
title('h as a function of t')
xlabel('t');
ylabel('h')




%%
%2.c
x0 = [1 0 0];
r1 = 0.04;
r2 = 10^4;
r3 = 3*10^7;

options = odeset('RelTol', 10^-6);
f = @(t, x)[-r1*x(1)+(r2*x(2)*x(3)); (r1*x(1))- (r2*x(2)*x(3))-(r3*(x(2)^2)); r3*(x(2)^2)];
[t, x] = ode23s(f, [0,1000], x0', options);
%867 utan
%10-3 : 31
%10-4 : 38
%10-5 : 49
%10-6 : 62

h=0;
for i = 1 : length(t)-1
    h(i) = abs(t(i) - t(i+1));
end

figure(12)
h = [h 0];
plot(t, h)

%%
%Assignment 3.a
r = 2;

func = @(x) [1-(r^2)*((x(1)^2-(x(2)^2))/(((x(1)^2)+(x(2)^2))^2)); (-2*x(1)*x(2)*(r^2))/(((x(1)^2)+(x(2)^2))^2)];

x0 = [-4 0.2];
start = cell(1,4);

Ui = cell(1,4);
Ti = cell(1,4);


yval = [0.2 0.6 1 1.6];
xval = -4;
for i = 1:4
    start{i} = [xval yval(i)];
end

for i = 1:4
    tspan = 0:.5:10;
    [t, u] = RK(func, start{i}, 2, 100, 10);
    Ui{i} = u';
    Ti{i} = t;
end
figure(12)
plot(Ui{1}(:,1), Ui{1}(:,2))
hold on;
plot(Ui{2}(:,1), Ui{2}(:,2))
hold on;
plot(Ui{3}(:,1), Ui{3}(:,2))
hold on;
plot(Ui{4}(:,1), Ui{4}(:,2))
title('Particle flow past a cylinder')
xlabel('x');
ylabel('y')
legend('1','2','3', '4')
axis equal
    

%%
%Assignment 3.b
H = 2;
v0 = 20;
alpha = [30, 45, 60];
k = 0.065;

f = @(y) [y(2); -k*y(2)*sqrt((y(2)^2)+(y(4)^2)); y(4); -9.81-k*y(4)*sqrt((y(2)^2)+(y(4)^2))];
Ulast = cell(1,3);
U2last = cell(1,3);
for i = 1:3
    u0 = [0, v0*cos(((alpha(i)*pi)/180)), H, v0*sin(((alpha(i)*pi)/180))];
    [t, u] = RK(f, u0, 4, 100, 10);
    u = u';
    Ulast{i} = u(:,3);
    U2last{i} = u(:,1);
    
end

figure(13)
plot(U2last{1}, Ulast{1},'-*','MarkerIndices',1:5:length(U2last{1}))
hold on
plot(U2last{2}, Ulast{2},'-*','MarkerIndices',1:5:length(U2last{2}))
hold on
plot(U2last{3}, Ulast{3},'-*','MarkerIndices',1:5:length(U2last{3}))
legend('30 degrees', '45 degrees', '60 degrees')
ylim([0 15])
title('motion of particle with different degrees')
xlabel('length');
ylabel('height')
%%

function [t, y] = RK(f, y0, inVals,N, len)


y=zeros(inVals,N+1);
y(:,1)=y0;
h=len/N;
t=linspace(0,len,N+1)';

for i = 1:N
    k1 = f(y(:,i));
    k2 = f(y(:,i)+h*k1);
    k3= f(y(:,i)+(h*k1/4)+(h*k2/4));
    y(:,i+1) = y(:,i) + ((h/6)*(k1+k2+4*k3));
end
end
