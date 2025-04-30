% DGSM_SIR
% ORIGINAL AUTHOR: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: February, 2025
%
% Uses the trajectory algorithm to identify Morris' indices
%%
clear; clc; %close all;
%% Initial conditions for the SIR model
S0 = 900; I0 = 100; R0 = 0;
X0 = [S0; I0; R0];
N = sum(X0);

% Time variables
t_final = 50;
dt = 0.05;
t_data = 0:dt:t_final;
ode_options = odeset('RelTol',1e-6);

% We will consider two quantities of interest: the integral of the infected
% state and the steady-state equilibrium value of the infected state
QoI_time = length(t_data);
Y0 = [S0; I0; R0];

% Identify upper and lower bounds for [gamma, k, r, mu]
parameter_names = {'$\gamma$','$k$','$r$','$\mu$'};
UB = [0.3 1.0 1.0 0.5];
LB = [0.1 0.1 0.1 0.1];

par_all = 0.5.*(UB+LB);
%%

%%
par_nom = par_all;
pars = par_all;
num_par = length(UB);


%%
M = 1000;           % Number of samples we want
p = num_par;   % Number of parameters
delta = 1e-6;  % Step Size
upper = UB(:)';
lower = LB(:)';
d = zeros(M,p); % Store the elementary effects
%% 
F_storage = zeros(p+1,M);
% F_storage = zeros(M,1);

mean_value = 0;
mean_sq_value = 0;
% lhs_sample = lhsdesign(M,p);
% qcurr = LB+(UB-LB).*lhs_sample;
qcurr = zeros(M,p);
for j=1:4
    qcurr(:,j) = unifrnd(LB(j)+delta,UB(j)-delta,1,M);
end
for i=1:M
    
    [~,fbase] = ode45(@SIR_model,t_data,Y0,ode_options,qcurr(i,:),N);
    F_storage(1,i) = fbase(end,3); % Store the solution in case we want to look at time series
    mean_value = mean_value+fbase(end,3);
    mean_sq_value = mean_sq_value+fbase(end,3).^2;
    for j=1:p
        par_plus = qcurr(i,:);
        par_plus(j)=par_plus(j)+delta;
        [~,fstep] = ode45(@SIR_model,t_data,Y0,ode_options,par_plus,N);
        % Calculate the DGSM
        d(i,j) =  (fstep(end,3) - fbase(end,3))./delta;

        F_storage(j,i) = fstep(end,3);

         mean_value = mean_value+fstep(end,3);
        mean_sq_value = mean_sq_value+fstep(end,3).^2;
    end
end

% Calculate the variance of all the outputs
mean_value = mean_value./(M*(p+1));
mean_sq_value = mean_sq_value./(M*(p+1));

mu = sum(d,1)./M;
mu_star= sum(abs(d),1)./M;
v = sum(d.^2,1)./M;
% E[X^2] - E[X]^2
var_out = mean_sq_value - mean_value.^2;
%%
figure(10);clf;hold on;
for i=1:num_par
    plot(mu_star(i),v(i),'k.','MarkerSize',20)
    text(mu_star(i)+20,v(i),parameter_names{i},'FontSize',16,'Interpreter','latex');
end

% Parameters are uniform, so poincare constant is (b-a).^2./pi.^2
poincare_const = (UB-LB).^2 ./ (pi.^2);
figure(20);clf;
bar(v);
xticklabels(parameter_names);

figure(40);clf;
bar(poincare_const.*v./var_out);
xticklabels(parameter_names);
%%

