clc
clear all

%% 2a



%[X,FLAG] = PCG(A,B,...) also returns a convergence FLAG:
%    0 PCG converged to the desired tolerance TOL within MAXIT iterations
%    1 PCG iterated MAXIT times but did not converge.
%    2 preconditioner M was ill-conditioned.
%    3 PCG stagnated (two consecutive iterates were the same).
%    4 one of the scalar quantities calculated during PCG became too
%      small or too large to continue computing.

%   [X,FLAG,RELRES] = PCG(A,B,...) also returns the relative residual
%   NORM(B-A*X)/NORM(B). If B is zero, then RELRES is set to 0, unless
%   the 'null' option is set, where the relative residual is defined as
%   NORM(A*X)/NORM(X). If FLAG is 0, then RELRES <= TOL.

%   [X,FLAG,RELRES,ITER,RESVEC] = PCG(A,B,...) also returns a vector of the
%   estimated residual norms at each iteration including NORM(B-A*X0) if
%   the 'null' option in absent. With the 'null' option, the residuals and
%   the relative residuals are defined to be the same.

load('cooling_flange.mat');
spy(A);

b = rand(35296,1);
tol = 10^-4;
MAXIT = 1000;

tic

[X1, FLAG, RELRES, ITER, RESVEC] = pcg(A,b,tol, MAXIT);

timeElapsed1 = toc

semilogy(0:ITER,RESVEC/norm(b), 'k.', 'MarkerSize', 5);
xlabel('iteration number, k')
ylabel('RESVEC')

ITER;
RELRES;

tic
x2 = A\b;

timeElapsed2 = toc

%% 2b
load('cooling_flange.mat');

b = rand(35296,1);
tol = 10^-4;
MAXIT = 1000;

tic
M = diag(diag(A)); 

[X1, FLAG, RELRES1, ITER1, RESVEC1] = pcg(A,b,tol, MAXIT, M); 
ITER1
RELRES1
timeelapsed1 = toc
semilogy(0:ITER1,RESVEC1/norm(b), 'b.', 'MarkerSize', 5);
xlabel('iteration number, k')
ylabel('RESVEC')
hold on

tic
L = ichol(A);

[X2, FLAG, RELRES2, ITER2, RESVEC2] = pcg(A,b,tol, MAXIT,L,L');
ITER2
RELRES2
timeelapsed2 = toc

semilogy(0:ITER2,RESVEC2/norm(b), 'r.', 'MarkerSize', 5);
legend( 'precon. conjugate gradient','incomplete Cholesky factorization')


%% 2c
%load('convdiff.mat')
load('cooling_flange.mat');
b = rand(55096,1);
tol = 10^-1; %%kan ej visa att den ej convergerar
MAXIT = 100;
% 
% [X1, FLAG] = pcg(A,b, 10^-12, 1000);
% FLAG;

% [X2, FLAG, RELRES, ITER, RESVEC] = gmres(A,b,1e-8, 100);
% tic
tic

A = sparse(A);
% p = symrcm(A);

setup.type = 'crout';
setup.milu = 'row';
setup.droptol = 10^-6;
% [L,U] = ilu(A,struct('type','ilutp','droptol',1e-6));
% [L,U] = ilu(A(p,p), setup);
L = sparse(L);
U = sparse(U);
% % spy(L)
% spy(U)
% 
[X2, FLAG2, RELRES, ITER, RESVEC] = gmres(A,b,tol,MAXIT,L,U);
FLAG2
% timeelapsed3 = toc
%%

%2ca

load('convdiff.mat');

b = rand(55096,1);
tol = 10^-4;
MAXIT = 1000;

tic

[X1, FLAG, RELRES, ITER, RESVEC] = gmres(A,b,[],tol, MAXIT);

timeElapsed1 = toc

semilogy(0:ITER(2),RESVEC/norm(b), 'k.', 'MarkerSize', 5);
title('GMRES tol 1e⁻4')
xlabel('iteration number, k')
ylabel('RESVEC')
tic
x = A\b;
timer = toc

ITER;RESVEC;
RELRES;

%%
%2cb
load('convdiff.mat');
b = rand(55096,1);
tol = 10^-4; %%kan ej visa att den ej convergerar
MAXIT = 1000;

tic
[L,U] = ilu(A);
[X2, FLAG2, RELRES, ITER, RESVEC] = gmres(A,b,[],tol,MAXIT,L,U);
timer = toc
FLAG2
semilogy(0:ITER(2),RESVEC/norm(U\(L\b)), 'k.', 'MarkerSize', 5);
title('GMRES tol 1e⁻4 with ilu')
xlabel('iteration number, k')
ylabel('RESVEC')
hold on
M = diag(diag(A));
tic
[X1, FLAG, RELRES1, ITER1, RESVEC1] = gmres(A,b,[],tol, MAXIT, M); 
timer = toc 
semilogy(0:ITER1(2),RESVEC1/norm((M\b)), 'b.', 'MarkerSize', 5);
legend('incomplete LU factorization', 'Precon using diag(A) = M')

%% pcg 
load('convdiff.mat');

b = rand(55096,1);
tol = 10^-4;
MAXIT = 300;

tic

[X1, FLAG, RELRES, ITER, RESVEC] = pcg(A,b,tol, MAXIT);

timeElapsed1 = toc
FLAG

ITER;RESVEC;
RELRES;
