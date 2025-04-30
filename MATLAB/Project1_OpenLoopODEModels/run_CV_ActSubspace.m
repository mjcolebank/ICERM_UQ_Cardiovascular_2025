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
%% Run active subspaces depending on the quantity of interest
f = @(q) act_subspace_f(q,IC,tspace,1,[]);
param_ids = 1:11; % Exclude the cardiac cycle length;
UB = param.*1.5;
LB = param.*0.5;
M = 100;
CS_flag = 1; % Use complex step, otherwise use centered finite difference
parallel_flag=1;
[C_f,D,W,act_scores] = ActiveSubspace_Algorithm(f,UB,LB,M,CS_flag,param_ids,param_base,parallel_flag);
[mu,mu_star,sigma] = Morris_Screening_Algorithm(f,UB,LB,M,param_ids,param_base)
%%
figure(10);clf;
subplot(1,2,1);
bar(D./max(D));
xticks(1:length(param_ids));
xticklabels(param_names{param_ids});
title('Eigenvalues'); grid on;
set(gca,'FontSize',20);

subplot(1,2,2);
bar(act_scores./max(act_scores))
xticks(1:num_param)
title('Activity Scores'); grid on;
set(gca,'FontSize',20);
%%
function out = act_subspace_f(q,IC,tspace,outflag,data)
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
    out = sum([res1 res1 res3].^2);
end

end
