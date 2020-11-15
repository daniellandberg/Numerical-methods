% Exercise 4 Huperbolic equations
% Upwind a > 0 stable if |a|lamdba =< 1
tmax = 6;
xmin = 0;
xmax = 4.5;
tau = 2.2;
N = 100;
dt =0.02;
dx = (xmax-xmin)/N;
lambda = dt/dx;
tmin = 0;
T_in = tmin : dt : tmax- dt;
a = 1;
%R=@(t) sin(2*pi*t/tau);
R=@(t) square(2*pi*t/tau);
gsq = square(2*pi*T_in/tau);
%gsin = sin(2*pi*T_in/tau);
u = zeros((round(tmax/dt)),(xmax/dx));
X_in=(xmin : dt : xmax-dt);
% BC
u(:,1) = gsq;
%u(:,1) = gsin;
tempVec = zeros(round(tmax/dt),1);
for n = 1 : round(tmax/dt) -1 %Time loop:
    for j = 2 : xmax/dx  % Space loops, Starting from 2, IC at 1
        u(n+1,j) = u(n,j) - a*lambda*(u(n,j)-u(n,j-1)); % upwind
    end
    A = circshift(tempVec,1);
    A(1) = R(T_in(n+1));
    tempVec = A;
    test = A(1:length(X_in));
    %test(n) = R(T_in(n));
    %plot((xmin : dx : xmax-dx),test(n))
    %plot((xmin : dx : xmax-dx),test, 'g')
    
    %plot((xmin : dx : xmax-dx),u(n,:),'r')
%     axis([0 xmax -1 1])
%     shg
%     pause(dt/10)
    
end
Upwind = u(end,:);
plot(X_in, test,'b')
hold on
plot((xmin : dx : xmax-dx), u(end,:),'r')
xlabel('x')
title('Plot at time = 6')
legend('exact solution', 'Upwind')
axis([0 xmax -1.1 1.1])
%%
% Lax?Wendroff, stable for |a|lambda =<1
tmin = 0;
tmax = 6;
xmin = 0;
xmax = 4.5;
tau = 2.2;
N = 100;
dt = 0.02; %0.03
dx = (xmax-xmin)/N;
lambda = dt/dx;
T_in = tmin : dt : tmax- dt;
a = 1;
X_in=(xmin : dt : xmax-dt);
%gsq = square(2*pi*T_in/tau);
gsin = sin(2*pi*T_in/tau);
R=@(t) sin(2*pi*t/tau);
%R=@(t) square(2*pi*t/tau);

u = zeros((round(tmax/dt)),(xmax/dx));
test = zeros((round(tmax/dt)),1);
% BC
%u(:,1) = gsq;
u(:,1) = gsin;
tempVec = zeros(round(tmax/dt),1);

for n = 1 : round(tmax/dt)-1 %Time loop:
    for j = 2 : xmax/dx  % Space loops, Starting from 2, IC at 1
        if j == xmax/dx %Extrapolate for j+1
            polate = 2*u(n, j-1) - u(n, j-2);
            u(n+1,j) = u(n,j) - a*lambda/2*(polate-u(n,j-1))+((a^2*lambda^2)/2)*(polate-2*u(n,j)+u(n,j-1));
        else
            u(n+1,j) = u(n,j) - a*lambda/2*(u(n,j+1)-u(n,j-1))+((a^2*lambda^2)/2)*(u(n,j+1)-2*u(n,j)+u(n,j-1)); % upwind
        end
        
    end
    A = circshift(tempVec,1);
    A(1) = R(T_in(n+1));
    tempVec = A;
    test = A(1:length(X_in));
%     plot((xmin : dx : xmax-dx),u(n,:),'r')
%     axis([0 xmax -1 1])
%     shg
%     pause(dt/1000)
%     
end
Wendroff=u(end,:);
plot(X_in, test,'b')
hold on
plot((xmin : dx : xmax-dx), u(end,:),'r')
xlabel('x')
title('Plot at time = 6')
legend( 'exact solution','Lax wendroff')
axis([0 xmax -1.1 1.1])

%%
% Lax?Friedrich, stable for |a|lambda =<1

tmax = 6;
xmin = 0;
xmax = 4.5;
tau = 2.2;
N = 100;
dt = 0.02;
dx = (xmax-xmin)/N;
lambda = dt/dx;
tmin = 0;
T_in = tmin : dt : tmax- dt;
X_in=(xmin : dt : xmax-dt);
tmax = 6;
a = 1;
R=@(t) sin(2*pi*t/tau);
%R=@(t) square(2*pi*t/tau);
%gsq = square(2*pi*T_in/tau);
gsin = sin(2*pi*T_in/tau);
u = zeros((round(tmax/dt)),(xmax/dx));
% BC
%u(:,1) = gsq;
u(:,1) = gsin;
tempVec = zeros(round(tmax/dt),1);

for n = 1 : round(tmax/dt)-1 %Time loop:
    for j = 2 : xmax/dx  % Space loops, Starting from 2, IC at 1
        if j == xmax/dx %Extrapolate for j+1
            polate = 2*u(n, j-1) - u(n, j-2);
            u(n+1,j) = 0.5*(polate+u(n,j-1)) - a*lambda/2*(polate-u(n,j-1));;
        else
            u(n+1,j) = 0.5*(u(n,j+1)+u(n,j-1)) - a*lambda/2*(u(n,j+1)-u(n,j-1)); % upwind
        end
    end
    A = circshift(tempVec,1);
    A(1) = R(T_in(n+1));
    tempVec = A;
    test = A(1:length(X_in));
%     plot((xmin : dx : xmax-dx),u(n,:),'r')
%     axis([0 xmax -1 1])
%     shg
%     pause(dt/10)
    
end
Friedrich=u(end,:);
plot(X_in, test,'b')
hold on
plot((xmin : dx : xmax-dx), u(end,:),'r')
xlabel('x')
title('Plot at time = 6')
legend( 'exact solution','Lax Friedrich')
axis([0 xmax -1.1 1.1])
%%
% Plotting the three methods above on top of each other
plot((xmin : dx : xmax-dx), Upwind,'g')
hold on
plot((xmin : dx : xmax-dx), Friedrich,'b')
hold on 
plot((xmin : dx : xmax-dx), Wendroff,'r')
hold on
plot(X_in, test,'--')
xlabel('x')
title('Plot at time = 6')
legend( 'Upwind','Lax Friedrich','Lax wendroff','Exact Solution')
axis([0 xmax -1.4 1.4])
%%
% part 2 upwind
tmin = 0;
tmax = 6;
xmin = 0;
xmax = 5;
alpha = 0.5;
Tcool = 50;
Thot = 200;
a = 1;
N = 100;
dt = 0.03; %0.03 gives exact = wandroff??
dx = (xmax-xmin)/N;
lambda = dt/dx;
% T_in = tmin : dt : tmax;
T_in = tmin : dt : tmax- dt;
BC1 = Tcool +(Thot - Tcool)*sin(4*pi*T_in);
BC2 = Thot;
BC3 = Thot + Tcool*sin(5*pi*(T_in-1));
T = zeros((round(tmax/dt)),(xmax/dx)); % IC

% Setting up BC
for i = 1 : length(T_in)
    if T_in(i) <=0.125
        T(i,1) = BC1(i);
    elseif (T_in(i) > 0.125) && (T_in(i) <= 1)
        T(i,1) = BC2;
    else
        T(i,1) = BC3(i);
    end
end
T(1,:) = ones(length(T(1,:)),1)*Tcool;
for n = 1 : round(tmax/dt)-1 %Time loop:
    for j = 2 : xmax/dx  % Space loops, Starting from 2, IC at 1
        T(n+1,j) = T(n,j) - a*lambda*(T(n,j)-T(n,j-1))-dt*(alpha*(T(n,j)-Tcool)); 
    end

%     mesh(T_in,(xmin : dx : xmax-dx), T')
%     shg
%     pause(dt/10)
    
end
T1 = T;
mesh(T_in,(xmin : dx : xmax-dx), T')
xlabel('t');
ylabel('x')
zlabel('Temperature')
%plot((tmin : dt : tmax), u(:,100))

%%
% part 2 Laxen
tmin = 0;
tmax = 6;
xmin = 0;
xmax = 5;
alpha = 0.5;
Tcool = 50;
Thot = 200;
a = 1;
N = 100;
dt = 0.03;
v = 1;
dx = (xmax-xmin)/N;
lambda = dt/dx;
T_in = tmin : dt : tmax- dt;
BC1 = Tcool +(Thot - Tcool)*sin(4*pi*T_in);
BC2 = Thot;
BC3 = Thot + Tcool*sin(5*pi*(T_in-1));
T = zeros((round(tmax/dt)),(xmax/dx)); % IC

% Setting up BC
for i = 1 : length(T_in)
    if T_in(i) <=0.125
        T(i,1) = BC1(i);
    elseif (T_in(i) > 0.125) && (T_in(i) <= 1)
        T(i,1) = BC2;
    else
        T(i,1) = BC3(i);
    end
end
sig = v*dt/dx;
T(1,:) = ones(length(T(1,:)),1)*Tcool;
for n = 1 : round(tmax/dt)-1 %Time loop:
    for j = 2 : xmax/dx  % Space loops, Starting from 2, IC at 1
        if j == xmax/dx % Extrapolate for j+1
            polate = 2*T(n, j-1) - T(n, j-2);
            T(n+1,j) = T(n,j) -(sig*(1-alpha*dt)/2)*(polate-T(n,j-1))+(0.5*sig^2)*(polate-2*T(n,j)+T(n,j-1))-dt*(1-0.5*alpha*dt)*(alpha*(T(n,j)-Tcool));
        else
            T(n+1,j) = T(n,j) -(sig*(1-alpha*dt)/2)*(T(n,j+1)-T(n,j-1))+(0.5*sig^2)*(T(n,j+1)-2*T(n,j)+T(n,j-1))-dt*(1-0.5*alpha*dt)*(alpha*(T(n,j)-Tcool));
        end
    end

%     plot(T_in,(xmin : dx : xmax-dx),T(n,:))
%     axis([0 xmax -1 1])
%     shg
%     pause(dt/10)
    
end
mesh(T_in,(xmin : dx : xmax-dx), T')
xlabel('t');
ylabel('x')
zlabel('Temperature')
%plot((tmin : dt : tmax), u(:,100))
%%
% subplots a)
subplot(1,2,1);
mesh(T_in,(xmin : dx : xmax-dx), T1')
title('Upwind')
xlabel('t');
ylabel('x')
axis tight
subplot(1,2,2);
mesh(T_in,(xmin : dx : xmax-dx), T')
title('Lax?Wendroff')
xlabel('t');
ylabel('x')
axis tight
%%
% b)
% at t = 3 and t = 6
idx = (tmax/dt);

plot((xmin : dx : xmax-dx),T(idx,:))
hold on
plot((xmin : dx : xmax-dx),T1(idx,:)) 
legend('Lax?Wendroff', 'Upwind')
xlabel('x')
ylabel('Temperature')
title('Plot showing the Tempature at different x values at t = 6')