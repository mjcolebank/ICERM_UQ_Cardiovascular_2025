% Predicts continuous and bolus dose concentrations of cytokines and
% monocytes to a bolus and continous dose of LPS. This computer code
% available under MIT open license as described in the readme.
% Developers: Kristen Windoloski and Mette Olufsen (NCSU).
% Contact: Mette Olufsen (msolufse@ncsu.edu)
%
clear all; close all; clc;

global dataB dataC

% -------- Continuous infusion data and model solution

% Set figure size
plotpars.height = 600;
plotpars.width  = 1400;
plotpars.lw     = 3; 
plotpars.fs     = 14;
xlim = [0 12];
    
% Load mean data and set up figure
meanctsdataplots(1,plotpars); % Set up figure
load 'MeanCtsOptPars.mat';    % Load optimized parameter values

% Run model for mean data
dataC.pars  = optpars; % Pptimized parameter values

dataC.ic    = cts_model_ic;     % Assigns initial conditions
dataC.tspan = 0:0.01:12;        % Set time span
sol  = cts_model_solver(dataC); % Solves ODE model
t    = dataC.tspan;             % Extract time solution
y    = deval(sol,t);            % Evaluate solution at specified times

% Plot model predictions with data
[ctsR2,ctsCurve] = cts_model_plotsOPT(1,t,y,sol,xlim,'r',plotpars);


% -------- Bolus data and model solution

% Load mean data and set up figure
meanbolusdataplots(1,plotpars); % Set up figure
load 'MeanBolusOptPars.mat'     % Load optimized parameter values

% Run model for mean data
dataB.pars  = optpars; % Optimized parameter values

dataB.ic    = bolus_model_ic;     % Assigns initial conditions
dataB.tspan = 0:0.01:12;          % Set time span
sol  = bolus_model_solver(dataB); % Solves ODE model
t    = dataB.tspan;               % Extract time solution
y    = deval(sol,t);              % Evaluate solution at specified times

%t_sol = t;
%y_sol = y;

% Plot model predictions with data
[bolusR2,bolusCurve] = bolus_model_plotsOPT(1,t,y,sol,xlim,'k',plotpars);

% Add R2 values
figure(1); 
subplot(2,4,5); 
legend('',['R^2 = ' num2str(ctsR2(1),'%.2f')],'',['R^2 = ' num2str(bolusR2(1),'%.2f')]);
subplot(2,4,6);
legend('',['R^2 = ' num2str(ctsR2(2),'%.2f')],'',['R^2 = ' num2str(bolusR2(2),'%.2f')]);
subplot(2,4,7);
legend('',['R^2 = ' num2str(ctsR2(3),'%.2f')],'',['R^2 = ' num2str(bolusR2(3),'%.2f')]);
subplot(2,4,8);
legend('',['R^2 = ' num2str(ctsR2(4),'%.2f')],'',['R^2 = ' num2str(bolusR2(4),'%.2f')]);
% exportgraphics(gcf,'optmeanmodel.eps');

