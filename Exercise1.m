clc
clear all
legend('10','20','40', '80', '160', '320')
N = [10,20,40,80,160,320]; 
T = cell(1,length(N));
U = cell(1,length(N));

f=@(y) [y(2); (1-y(1)^2)*y(2)-y(1)];

for i = 1 : length(N)

    h = 1/N(i); %Step Length 
    [t, u] = RK(f, [1 0], 2, N(i));
    T{i} = t;
    U{i} = u;
    figure(1);
    plot(t,u(1,:)');
    hold on;
    figure(1)
    xlabel('time');
    title('ex 1')
    hold on
    
end
legend('10','20','40', '80', '160', '320')
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
figure(2);
loglog(hLista,en);
xlabel('h')
ylabel('en')


%% 

%ex21
r1 = 0.04;
r2 = 10^4;
r3 = 3*10^7;
f = @(x)[-r1*x(1)+(r2*x(2)*x(3)); (r1*x(1))- (r2*x(2)*x(3))-(r3*(x(2)^2)); r3*(x(2)^2)];
N = [125, 250, 500, 1000, 2000];
Steg = [500, 510, 520, 530, 540, 550, 560, 570, 580, 580, 600]; %530 är det minsta N som ger stabilitet
x0 = [1 0 0];

for i = 1 : length(N)
   [t, u] = RK(f, x0, 3,N(i));
  
    T2{i} = t;
    U2{i} = u;
    
end
u  = U2{5};
t = T2{5};

figure(4);
loglog(t,u(1,:)');
hold on;

figure(5)
loglog(t,u(2,:)');
hold on;

figure(6)
loglog(t,u(3,:)');

en = U2{5}(:,end)-U2{4}(:,end);

%%
%%22
x0 = [1 0 0];
r1 = 0.04;
r2 = 10^4;
r3 = 3*10^7;
options = odeset('RelTol', 10^-6);
f = @(t, x)[-r1*x(1)+(r2*x(2)*x(3)); (r1*x(1))- (r2*x(2)*x(3))-(r3*(x(2)^2)); r3*(x(2)^2)];
[t, x] = ode23(f, [0,1], x0', options);
%figure(10)
%plot(t,x(:,2));

%867 utan
%10-3 : 867
%10-4 : 868
%10-5 : 869
%10-6 : 869

h=0;
for i = 1 : length(t)-1
    h(i) = abs(t(i) - t(i+1));
end

figure(11)
h = [h 0];
plot(t, h)




%%
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

