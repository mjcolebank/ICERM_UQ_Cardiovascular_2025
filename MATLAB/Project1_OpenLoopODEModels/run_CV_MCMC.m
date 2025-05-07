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
load CV_test_data.mat data param_true
%% Run active subspaces depending on the quantity of interest
param_ids = 1:11; % Exclude the cardiac cycle length;
theta0 = param_base(param_ids);
UB = theta0.*1.5;
LB = theta0.*0.5;
f = @(q) MCMC_f(q,IC,tspace,1,data,param_ids,param_base);


prior_F = @(dummy) 1./prod(UB-LB); % The prior PDF is constant throughout the domain (equiprobable outcomes
M = 1000;
k0 = 100;
num_par = length(param_ids);
covar = eye(num_par,num_par).*0.01.*abs(theta0);
[chain,s2chain] = AdaptiveMetropolis(f,data,prior_F,theta0,k0,M,covar,UB,LB);
%%
figure(1);clf;
if num_par<5
    for i=1:num_par
        subplot(2,2,i);
        plot(chain(i,:));
        yline(param_true(param_ids(i)),'--k');
    end
elseif num_par<10
    for i=1:num_par
        subplot(3,3,i);
        plot(chain(i,:));
        yline(param_true(param_ids(i)),'--k');
    end
elseif num_par<13
    for i=1:num_par
        subplot(3,4,i);
        plot(chain(i,:));
        yline(param_true(param_ids(i)),'--k');
    end
end

figure(2);clf;
plot(s2chain);
%%

%%
function out = MCMC_f(q,IC,tspace,outflag,data,param_ids,param_base)
param = param_base;
param(param_ids) = q;
[output,~] = call_CV_model(IC,param,tspace);

if outflag==1 % LV data
    res = [(max(output.plv) - max(data.plv))./max(data.plv) ...
           (min(output.plv) - min(data.plv))./min(data.plv)];
    out = sum(res.^2);
elseif outflag==2 %LV and Ao pressure
    res = [(max(output.plv) - max(data.plv))./max(data.plv) ...
           (min(output.plv) - min(data.plv))./min(data.plv) ...
           (max(output.pao) - max(data.pao))./max(data.pao) ...
           (min(output.pao) - min(data.pao))./min(data.pao)];
    out = sum(res.^2);
elseif outflag==3 %LV + Ao time series + Ao static
    res = [(max(output.plv) - max(data.plv))./max(data.plv) ...
           (min(output.plv) - min(data.plv))./min(data.plv) ...
           (max(output.pao) - max(data.pao))./max(data.pao) ...
           (min(output.pao) - min(data.pao))./min(data.pao) ...
           (output.pao - data.pao)./mean(data.pao)./sqrt(50)];
    out = sum(res.^2);
elseif outflag==4 %residual function (assume pressure in the Aorta, aortic and mitral flow
    res1 = (data.pao - output.pao)./max(data.pao);
    res2 = (data.qav - output.qav)./max(data.qav);
    res3 = (data.qmv - output.qmv)./max(data.qmv);
    out = sum([res1 res2 res3].^2);
end

end
