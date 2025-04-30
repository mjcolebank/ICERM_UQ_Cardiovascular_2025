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

outflag = 1;
%%
[IC, param] = get_CV_parameters();
T = param(12);

% Starting and ending time
tstart = 0;
tend   = 30.*T; % 30 cycles
dt = 1e-3;
tspace = tstart:dt:tend;
param_base = param;

%% Generate synthetic data from the model
[output_true,~] = call_CV_model(IC,param,tspace);
nt = length(output_true.plv);
data.plv = output_true.plv+normrnd(0,1,1,nt);
data.pao = output_true.pao+normrnd(0,1,1,nt);
data.qmv = output_true.qmv+normrnd(0,8,1,nt);
data.qav = output_true.qav+normrnd(0,8,1,nt);

%% Try to infer the model parameters using OLS
f = @(q) get_CV_residual(q,IC,tspace,outflag,data);
res = f(param_base);
param_ids = 1:11; % Exclude the cardiac cycle length;
UB = param.*1.5;
LB = param.*0.5;
M = 100;
CS_flag = 1; % Use complex step, otherwise use centered finite difference
parallel_flag=1;
[mu,mu_star,sigma] = Morris_Screening_Algorithm(f,UB,LB,M,param_ids,param_base);
%%
figure(1);clf;hold on;
for i=1:length(param_ids)
    plot(mu_star(i),sigma(i),'k.','MarkerSize',20)
    text(mu_star(i).*1.01,sigma(i),param_names{i},'FontSize',16,'Interpreter','latex');
end

rank = sqrt(mu_star.^2 + sigma.^2);
figure(2); clf;
bar(rank);
xticklabels(param_names);
%%
function res = get_CV_residual(q,IC,tspace,outflag,data)
[output,~] = call_CV_model(IC,q,tspace);

if outflag == 1 %Static pressure only
    res = [(max(data.plv)-max(output.plv))./max(data.plv) ...
           (min(data.plv)-min(output.plv))./min(data.plv) ...
           (max(data.pao)-max(output.pao))./max(data.pao) ...
           (min(data.pao)-min(output.pao))./min(data.pao)];
elseif outflag==2 %Full aortic profile (assume pressure in the Aorta)
    res = (data.pao - output.pao)./data.pao;
elseif outflag==3 %Static + dynamic pressure measurements only
    res = [(max(data.plv)-max(output.plv))./max(data.plv) ...
           (min(data.plv)-min(output.plv))./min(data.plv) ...
            (data.pao - output.pao)./data.pao];
elseif outflag==4 %mixed continuous data
    res = [(data.pao - output.pao)./data.pao ...
           (data.qmv - output.qmv)./max(data.qmv) ...
           (data.qav - output.qav)./max(data.qav)];
elseif outflag==5 %mixed static + dynamic data
    res = [(data.pao - output.pao)./data.pao ...
        (max(data.qmv)-max(output.qmv))./max(data.qmv) ...
           (max(data.qao)-max(output.qao))./max(data.qao)];
end

end
%%
function out = get_CV_output(q,IC,tspace,outflag)
param = q;
[output,~] = call_CV_model(IC,param,tspace);

if outflag == 1 %Static pressure only
    out = [max(output.plv) min(output.plv) max(output.pao) min(output.pao)];
elseif outflag==2 %Full aortic profile (assume pressure in the Aorta)
    out = output.pao;
elseif outflag==3 %Static + dynamic pressure measurements only
    out = [max(output.plv) min(output.plv) output.pao];
elseif outflag==4 %mixed continuous data
    out = [output.pao output.qmv output.qav];
elseif outflag==5 %mixed static + dynamic data
    out = [output.pao max(output.qmv) max(output.qav)];
end

end
