%=============================================================*
%                                                             *
% run_1D.m                                                    *
% Version: 2.0 (created on 6 July 2021)                       *
% AUTHORS: M.J. Colebank, M.U. Qureshi,  M.S. Olufsen         *
%          Department of Mathematics                          *
%          North Carolina State University, Raleigh, USA      *
% DATE UPDATED: 2 May 2025.                                 *
%                                                             *
% DESCRIPTION: This script creates the inteface with the C++  *
% code for the 1D fluid dynamics model. The script creates    *
% and runs the executable by passing selected model           *
% parameters and plots the hemodynamic waveforms              *
%                                                             *
%=============================================================*



% Model usines a simple linear elastic wall model and a three element Windkessl outflow
% boundary condition. The system is driven by imposing an inflow profile.

clc; clear all; close all;
%% Run the lines below if you are compiling for the first time
% !chmod +x sor06
% !make clean
% !make
format shortg;
%% Define the parameters
% The stiffness is assumed to follow an exponential curve, i.e.
% Eh/r0 = k1*exp(k2*r0) + k3
% where E, h, and r0 are the Youngs Modulus, wall thickness, and reference
% radius, respectively. To set a single stiffness value for all the
% vessels, set k1=0;
k1 = 0;%1e+5;
k2 = -25;
k3 = 8e4;%5e+4;


%% Specify Windkessel boundary
Rp_1 = 1;
Rp_2 = 1;
Rd_1 = 1.5;
Rd_2 = 1.5;
CT_1  = 1.4;
CT_2  = 1.4;

%% Define scaling factors for resistance/compliance if needed
pars = [k1,k2,k3,Rp_1,Rp_2,Rd_1,Rd_2,CT_1,CT_2,1];
pars_str = mat2str(pars);
%% call the model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unix(sprintf('sor06.exe  %s',pars_str(2:end-1)));
% For MJ Colebank laptop to run
% f_file = fopen('run_code.txt','w');
% fprintf(f_file,'make clean\n');
% fprintf(f_file,'make\n');
% fprintf(f_file,'chmod +x sor06\n');
% fprintf(f_file,'./sor06.exe %s\n',pars_str(2:end-1));
% fclose(f_file);
% return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting Simulated data
% Plot inlet/outlet waveforms in the whole tree
% data = load('art_ALL.2d');
data = dlmread('output_1.2d');
[time,x,p,q,A,C] = gnuplot(data);
t = time(:,1)-time(1,1); % Time starts from 0
% Proximal predictions are 1->N
% Midpoint predictions are N+1->2N
% Distal   prdictions are 2N+1->3N (or end)
figure(1); clf; hold on;
% plot(linspace(0,T,length(Pdat)),Pdat,'k','LineWidth',3);
plot(t,p(:,1),'LineWidth',2);
set(gca, 'FontSize',30);grid on;
xlabel('Time (s)');
ylabel('Pressure (mmHg)');
hold off;


