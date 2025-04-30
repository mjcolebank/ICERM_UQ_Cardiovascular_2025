% run_cardiovascular.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
% Defines initial conditions, model parameters, and other components needed
% to solve the 2-component cardiovascular model.
close all; clear; clc;
%%
[IC, param] = get_CV_parameters();
T = param(12);

% Starting and ending time
tstart = 0;
tend   = 50.*T; % 30 cycles
dt = 1e-3;
tspace = tstart:dt:tend;


[output,tplot] = call_CV_model(IC,param,tspace);

%% Plotting
% Pressure
figure(1); clf; hold on;
plot(tplot,output.plv,':c','LineWidth',2);
plot(tplot,output.pao,':m','LineWidth',2);
yline([120 80],'--k')
ylabel('Pressure (mmHg)')
xlabel('Time (s)')
grid on;
set(gca,'FontSize',20)

%% PV loop
figure(2); clf; hold on;
plot(output.Vlv,output.plv,'r','LineWidth',3);
ylabel('Pressure (mmHg)')
xlabel('Volume (mL)')
grid on;
set(gca,'FontSize',20);

%% Flows
figure(3); clf; hold on;
plot(tplot,output.qart,'r','LineWidth',3);
plot(tplot,output.qav,'b','LineWidth',3);
plot(tplot,output.qmv,'c','LineWidth',3);
ylabel('Pressure (mmHg)')
xlabel('Volume (mL)')
grid on;
set(gca,'FontSize',20);
legend('q art','q av','q mv')

