%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;
% this programe is used to estimate the parameter of 
%                     0.5
% G(s) = -----------------    
%               s^2 + s + 1 
%                0.009335 z + 0.008732
% G(z) = ---------------------------------------------           ; T = 0.2
%                z^2 - 1.783 z + 0.8187
%% Plot option
% if Plot(1) = 0  --> no plot needed
% if Plot(1) = 1  -->  plot needed
Plot=[1 1:4];
%% input signal 
T=0.2;
t=0:T:1000;
mu_u=0.7; % input mean
sigma_u=6; % input variance
mu_e=1; % noise mean
sigma_e=1; % noise variance
u= random('normal', mu_u, sigma_u,length(t),1); % input
eps= random('normal', mu_e, sigma_e,length(t),1); % noise addition
%% system characteristics
A=[1 -0.4 0.5];
B=[1.2, 0.3];
%% noise characteristics
C=[1, 0.8 -0.1];
%% estimation parameters
na=2;
nb=1;
d=1;
nc=2;
%% output signal
y(1)=0;
yn(1)=0;
for k=1:length(u)-1
    [ y_output ] = OutputEstimation( A, B, 1, u(1:k+1), y(1:k), k+1 );
    y(k+1)=y_output;
    [ y_output ] = OutputEstimation( A, C, 1, eps(1:k+1), yn(1:k), k+1 );
    yn(k+1)=y_output;
end
y1=y+yn;
%% estimation
[ Guz, Gez ] = ModifiedExtendedRecursiveLeastSquareVaringExponentialForgetting( u, y1, na, nb, d, nc, T,  Plot)