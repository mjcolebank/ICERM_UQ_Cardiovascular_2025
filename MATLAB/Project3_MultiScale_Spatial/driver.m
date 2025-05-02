% Driver file for running the Triseg model
% Model based on Colebank & Chesler 2022; adapted for humans
% Author: MJ Colebank
% Last edited: 5/1/2025
%
% THIS CODE COMES WITH NO GUARANTEES
% Abbreviations: LA -> left atrium, LV -> left ventricle, RA -> right
% atrium, RV -> right ventricle, S -> septum, Am -> midwall area, Cm ->
% midwall curvature, ef -> wall strain
clear; clc; close all;
%% First, load in nominal parameters and initial conditons
% Note that parameters are stored in a structure like format
[IC,pars] = get_pars_GR; % Get the initial conditions and parameter values
%% Next, define the time vector
num_cycles = 30;                 % Number of cycles to run the model for
dt         = 1e-4;               % Time stepping for the output and plotting
Tend       = num_cycles.*pars.T; % The end time point for the ODE/DAE solver
tspace     = 0:dt:Tend;          % The time vector
pars.dt    = dt;                 % Append the time stepping

%% Now solve the model
% This includes a system of differential algebraic equations (DAEs) which
% need to provide a mass matrix for enforcing contraints
%% MASS MATRIX DAE APPROACH
M = eye(20);  % Define the mass matrix
M(9,9)   = 0; % DAE - Vsw
M(20,20) = 0; % DAE - ym

options=odeset('Mass',M,'RelTol',1e-8, 'AbsTol',1e-8);  %sets how accurate the ODE solver is

%% Timing statements and model solve using ODE15s
tode = [0 0];
r_conv = 1e6;
r_stop = 1e-3;
i=1; yold = [];
while r_conv>=r_stop % || i<5
    % for i=1:num_cycles
    tode   = [tode(2) i.*pars.T];
    ysol   = ode15s(@DE_model_mass,tode,IC,options,pars);
    tspace = linspace(ysol.x(1),ysol.x(end),100);
    ynew = deval(ysol,tspace);
    if ~isempty(yold)
        r_conv = sum((ynew(:)-yold(:)).^2);
    end
    IC = ysol.y(:,end-1);
    yold = ynew;
    i=i+1;
    %         figure(1); hold on; plot(ynew(4,:))
end

%% Obtain all relevant outputs (volumes, pressures, sarcomere dynamics)
toutput = tspace;
youtput = ynew;

% Get model outputs
outputs = get_model_results(toutput,youtput,pars);


%% Now plot outputs
p      = outputs.p./0.133322;  % KPa  --> mmHg  (see get_model_results for output order)
V      = outputs.V;       % mL   --> muL   (see get_model_results for output order)
q      = outputs.q;       % mL/s --> muL/s (see get_model_results for output order)
Ca2    = outputs.Sarc(1:5,:);  % Contractility (Gamma in the manuscript) (LA,LV,RA,RV,S)
SL     = outputs.Sarc(6:10,:); % Sarcomere contractile length (Lsc)      (LA,LV,RA,RV,S)
Am     = outputs.Am;           % Midwall area of cardiac chambers        (LA,LV,RA,RV,S)
Cm     = outputs.Cm;           % Curvature of cardiac chambers           (LA,LV,RA,RV,S)
ef     = outputs.ef;           % Strain of the cardiac chamber           (LA,LV,RA,RV,S)
stress = outputs.stress;       % Wall stress of the cardiac chamber      (LA,LV,RA,RV,S)

% Time vector for plotting purposes
t      = toutput - toutput(1);
toutput = toutput-toutput(1);

% Double the vectors to see two cycles
t = [t t+t(end)];
p = repmat(p,1,2);
V = repmat(V,1,2);
q = repmat(q,1,2);
SL = repmat(SL,1,2);
ef = repmat(ef,1,2);
stress = repmat(stress,1,2);
Cm = repmat(Cm,1,2);
Am = repmat(Am,1,2);
%%
figure(1); clf;
subplot(2,2,1); hold on;
plot(V(7,:),p(7,:),'LineWidth',3);
title('LA');
set(gca,'FontSize',20);
ylabel('Pressure (mmHg');
subplot(2,2,2); hold on;
plot(V(3,:),p(3,:),'LineWidth',3);
title('RA');
set(gca,'FontSize',20);
subplot(2,2,3); hold on;
plot(V(8,:),p(8,:),'LineWidth',3);
ylabel('Pressure (mmHg');
xlabel('Volume (ml)');
title('LV');
set(gca,'FontSize',20);
subplot(2,2,4); hold on;
plot(V(4,:),p(4,:),'LineWidth',3);
title('RV');
xlabel('Volume (ml)');
set(gca,'FontSize',20);
%%
figure(2); clf;
subplot(1,2,1); hold on;
plot(t,p([7 8 1],:),'LineWidth',3); legend('LA','LV','SA');
ylim([0 max(p(1,:))])
ylabel('Pressure (mmHg)'); xlabel('Time (s)');
set(gca,'FontSize',20); %axis tight;
subplot(1,2,2); hold on;
plot(t,p([3 4 5],:),'LineWidth',3); legend('RA','RV','PA');
ylim([0 max(p(5,:))])
ylabel('Pressure (mmHg)'); xlabel('Time (s)');
set(gca,'FontSize',20); %axis tight;



%% 
long_strain = (SL-SL(:,1))./SL(:,1);
green_strain = 0.5.*(SL.^2./SL(:,1).^2 - 1.0);
circ_strain = (Am-Am(:,1))./Am(:,1);


figure(3);clf;
subplot(1,3,1);
plot(t,long_strain([2 4 5],:)','LineWidth',3)
grid on; set(gca,'FontSize',20);
ylim([-0.35 0.05])

subplot(1,3,2);
plot(t,circ_strain([2 4 5],:)','LineWidth',3)
grid on; set(gca,'FontSize',20);
ylim([-0.35 0.05])

subplot(1,3,3);
plot(t,green_strain([2 4 5],:)','LineWidth',3)
grid on; set(gca,'FontSize',20);
ylim([-0.3 0.05])
legend('LV','RV','S')


