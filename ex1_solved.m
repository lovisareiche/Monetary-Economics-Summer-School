%% REDS-SOLDS solution method for dynamic systems
% REDS-SOLDS is a package of Matlab codes written to solve rational expectations 
% models numerically, and to analyze the solution thus obtained. 
% Together, REDS.M and SOLDS.M are a simplified version of a package of codes 
% written by Robert King and Mark Watson, implementing the algorithms described 
% in their paper “System Reduction and Solution Algorithms for Singular Linear 
% Difference Systems Under Rational Expectations” (mimeo, 1995).
%
% REDS-SOLDS takes as input a model written in the form: A Et yt+1 = B yt + C xt
% where {xt} is a martingale difference sequence and yt is ordered so that 
% variables that are predetermined appear last in a subvector kt. We denote 
% NY = dim(yt), NX = dim(xt), and NK = dim(kt). We input A, B, C, NY, NX, NK 
%
% The program REDS.M reduces the system, i.e., transforms it so that it contains 
% a non-singular subsystem that can be solved and turned into a solution of 
% the whole model. This whole solution operation is performed by SOLDS.M, 
% whose output are the matrices D, F, G, and H in:
% yt = D kt + F xt
% kt+1 = G kt + H xt

clear 
close all

% assign parameters values - input for matrix construction
v.beta = 0.99; % discount factor
v.sigma = 1; % elasticity of intertemporal substitution: individuals are indifferent between consuming a certain amount of goods or services today and consuming the same amount in the future when the interest rate changes
v.varphi = 1; %  Frisch elasticity of labor supply: measures how responsive an individual's labor supply is to changes in their real wage rate
v.alpha = 0.67; % Calvo parameter: α is a measure of price stickiness (Every period a firm can reset its price with probability 1 − α)
v.phi_pi = 1.5; % policy parameter on inflation
v.phi_y = 0.5/4; % policy parameter on GDP
v.lambda = ((1-v.alpha)*(1-v.alpha*v.beta))/v.alpha; % parameter in NKPC
v.rho_a = 0.9;  % persistence of tech shock
v.rho_nu = 0.5;  % persistence of monpol shock
T = 13;     %number of periods for irfs


[A, B, C, Indicator_Variables, NY, NX, NK ]= Matrix_solved(v);
% Matrix solved rewrites the system of equations
% check Davids notes page 370
reds
solds

% keep real component of coefficient matrices if imaginary component small enough
if max(abs(imag(D)))<10^(-10)
    D=real(D);
end
if max(abs(imag(F)))<10^(-10)
    F=real(F);
end

%% Monetary policy shock
irf_nu = Irf_modif(2,13, D,F,G,H)

%compute manually with irf_modif results
res = zeros(5,14)
shock = [0 1]'
res(4:5,1) = H*shock  %a/nu in last two lines: compute response to shock and then gradual decay
for i=2:14
    res(4:5,i) = G*res(4:5,i-1) % G gives evolution of the predetermined vars (a and nu)
end
res(1:5,:) = D*res(4:5,:) % D gives evolution of control vars 

%% Technology shock
irf_a = Irf_modif(1,13, D,F,G,H)

%compute manually with irf_modif results
res = zeros(5,14)
shock = [1 0]'
res(4:5,1) = H*shock  %a/nu in last two lines
for i=2:14
    res(4:5,i) = G*res(4:5,i-1)
end
res(1:5,:) = D*res(4:5,:)
