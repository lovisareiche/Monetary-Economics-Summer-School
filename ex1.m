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

v.beta = 0.99;
v.sigma = 1;
v.varphi = 1; 
v.alpha = 0.67; 
v.phi_pi = 1.5;
v.phi_y= 0.5/4;
v.lambda = ((1-v.alpha)*(1-v.alpha*v.beta))/v.alpha; 
v.rho_a = 0.9;
v.rho_nu = 0.5; 
T=13; 


[A, B, C, Indicator_Variables, NY, NX, NK ]= Matrix(v);
reds % will yield error message: need to fill in Matrix first!
solds

% keep real component of coefficient matrices if imaginary component small enough
if max(abs(imag(D)))<10^(-10)
    D=real(D);
end
if max(abs(imag(F)))<10^(-10)
    F=real(F);
end

%% Monetary policy shock
%use inbuilt function 
irf_nu = Irf_modif(2,13, D,F,G,H)

%compute manually: create your own irf here

% Plot irf_nu for output, and your own irfs, do they look the same?

%% Technology shock
%use inbuilt function
irf_a = Irf_modif(1,13, D,F,G,H)

%compute manually: create your own irf here

% Plot irf_a for output, and your own irfs, do they look the same?