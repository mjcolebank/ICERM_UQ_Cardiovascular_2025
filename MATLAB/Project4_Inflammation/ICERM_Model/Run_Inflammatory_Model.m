% Predicts continuous and bolus dose concentrations of cytokines and
% monocytes to a bolus and continous dose of LPS. This computer code
% available under MIT open license as described in the readme.
% Developers: Kristen Windoloski and Mette Olufsen (NCSU).
% Contact: Mette Olufsen (msolufse@ncsu.edu)
% 
% LAST EDITED BY: Mitchel Colebank (mjcolebank@sc.edu)
% Edit date: 5/5/2025
clear all; close all; clc;

global dataB dataC

output_names = {'Endotoxin','Resting Macrophage',...
   'Activated Marcophage', 'TNF-$\alpha$', 'IL-6',...
   'IL-8','IL-10'};

load 'MeanCtsOptPars.mat';    % Load optimized parameter values

% Run model for mean data
dataC.pars  = optpars; % Pptimized parameter values

dataC.ic    = cts_model_ic;     % Assigns initial conditions
dataC.tspan = 0:0.01:12;        % Set time span
sol  = cts_model_solver(dataC); % Solves ODE model
t    = dataC.tspan;             % Extract time solution
y    = deval(sol,t);            % Evaluate solution at specified times

%%
figure(1);clf;
for i=1:7
    subplot(2,4,i);
    plot(t,y(i,:),'r','LineWidth',3);
    grid on; set(gca,'FontSize',20);
    title(output_names{i},'Interpreter','latex');
end


    
    