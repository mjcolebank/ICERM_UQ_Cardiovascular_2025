% run_cardiovascular.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: January, 2025
%
% Defines initial conditions, model parameters, and other components needed
% to solve the 2-component cardiovascular model.
close all; clear; clc;
addpath('../Project0_UQMethods/')
param_names = {'P LA', 'P SysCap', 'Rmv', 'Rav', 'Rart', 'Cao', 'Emax', 'Emin', 'Vlvd', 'T peak', 'T relax', 'T'};
%%
[IC, param] = get_CV_parameters();
T = param(12);

% Starting and ending time
tstart = 0;
tend   = 30.*T; % 30 cycles
dt = 1e-3;
tspace = tstart:dt:tend;
param_base = param;
load CV_test_data.mat data
%% Run active subspaces depending on the quantity of interest
% f = @(q) DGSM_f(q,IC,tspace,1,[]);
f = @(q) DGSM_f(q,IC,tspace,4,data); % For last QoI, use cost function

param_ids = 1:11; % Exclude the cardiac cycle length;
UB = param.*1.5;
LB = param.*0.5;
M = 10;
CS_flag = 1; % Use complex step, otherwise use centered finite difference
parallel_flag=0;
[mu,mu_star,v] = DGSM(f,UB,LB,M,param_ids,param_base,parallel_flag,CS_flag);
%%
% You can log scale v and mu_star if only a few parameters dominate
v = log(v); mu_star = log(mu_star);
%%
figure(1);clf;hold on;
for i=1:length(param_ids)
    plot(mu_star(i),v(i),'k.','MarkerSize',20)
    text(mu_star(i).*1.1,v(i),param_names{i},'FontSize',16,'Interpreter','latex');
end
set(gca,'FontSize',20);
grid on;

figure(2);clf;
bar(v);
xticks(1:length(param_ids));
xticklabels(param_names);
set(gca,'FontSize',20);
grid on;
%%

%%
function out = DGSM_f(q,IC,tspace,outflag,data)
param = q;
[output,~] = call_CV_model(IC,param,tspace);

if outflag==1 % systolic LV pressure
    out = max(output.plv);
elseif outflag==2 %systolic Ao pressure
    out = max(output.pao);
elseif outflag==3 %diastolic Ao pressure
    out = min(output.pao);
elseif outflag==4 %residual function (assume pressure in the Aorta, aortic and mitral flow
    res1 = (data.pao - output.pao)./max(data.pao);
    res2 = (data.qav - output.qav)./max(data.qav);
    res3 = (data.qmv - output.qmv)./max(data.qmv);
    out = sum([res1 res2 res3].^2);
end

end


