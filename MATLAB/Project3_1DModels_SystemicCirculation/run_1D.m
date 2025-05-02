%=============================================================*
%                                                             
% This is a driver file that will run the systemic fluids model from c++ and
% fortran by passing parameter values needed by the model.
%
% Authors: A LaPole, MJ Colebank, MS Olufsen
%=============================================================*

%%
clear; clc; close all;
%% Ensure the make file is compiled
!make clean
!make
% make the file executable
! chmod +x sor06

%% Define the connectivity matrix

% Converging full network
writeconn  	 = [0 1 2 0
                1 3 4 0
                2 0 0 0
                3 5 6 0
                4 0 0 0
                5 7 8 0
                6 0 0 0 
                7 0 0 0
                8 0 0 0];
            
terminal 	 = [2 4 6 7 8];

% Write all this information to file.
dlmwrite('connectivity.txt',writeconn,'\t');
dlmwrite('terminal_vessels.txt',terminal,'\t');

%% FOR PARTICIPANTS TO EDIT
% Consider some random values for the parameters of interest
% within the appropriate range
f2   = -35;
f3   = 400000;
fs2  = f2;
fs3  = f3;
alpha   = 0.90;

% File ID (useful for running the C++/Fortran code in parallel)
ID = 1;

%% Put the parameters into a vector and run the model
par_nom = [f2,f3,fs2,fs3,alpha,ID];
param_str = mat2str(par_nom);

% Run the model
% NOTE: Windows users need 'sor06.exe', Mac/Linux users need ./sor06
% out = unix(sprintf('./sor06 %s',param_str(2:end-1)));
out = unix(sprintf('sor06.exe %s',param_str(2:end-1)));

if out == 1
    disp 'there is a model output'
else
    disp 'there is no model output'
end

%% Load all the model results
% NOTE: Results are stored such that the first 1:N entries are the PROXIMAL
% large vessel predictions, N+1:2N are the
% MIDPOINT predictions in the same vessels, and 2N+1:3N are the DISTAL
% predictions.

name = sprintf('output_%d.2d',ID);

data = load(name);

%% Extract data of interest

% p - pressure
% q - flow
[~,~,p,q,~,~] = gnuplot(data); % extract data

nv = 9; % total no of vessels

ntp = size(p,1); % no of time points in the flow or pressure time series

pressure_all = (p(:,nv+1:2*nv)); % middle prediction all vessels
flow_all = (q(:,nv+1:2*nv)); % middle prediction all vessels

% competition data consists of flow (entire time series) from vessels 1-7, 
% and pressure (min and max points only) from vessel 9, 
% so these are the predictions we save
pressure_inference  = [max(pressure_all(:,9)), min(pressure_all(:,9))];
flow_inference = flow_all(:,1:7);

%% Do some plotting
% Time vector
t = linspace(0,0.7,ntp);

% Plot midpoint in the vessels
figure(1);clf(1);
plot(t,p,'LineWidth',3);
ylabel('Pressure (mmHg)');
xlabel('Time (s)');
grid on; set(gca,'FontSize',20);

figure(2);clf(2)
plot(t,q,'LineWidth',3);
ylabel('Flow (mL/s)')
xlabel('Time (s)');
grid on; set(gca,'FontSize',20);