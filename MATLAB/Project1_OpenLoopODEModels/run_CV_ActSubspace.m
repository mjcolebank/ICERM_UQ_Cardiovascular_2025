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
load('CV_test_data.mat','data')
% Starting and ending time
tstart = 0;
tend   = 20.*T; % 30 cycles
dt = 1e-3;
tspace = tstart:dt:tend;
param_base = log(param);
%% Run active subspaces depending on the quantity of interest
f = @(q) act_subspace_f(q,IC,tspace,4,data);
param_ids = [1:11];%1:11; % Exclude the cardiac cycle length;
UB = param_base.*1.2;
LB = param_base.*0.8;
% LB(4) = param(4).*0.95;
% UB(4) = param(4).*1.05;
M = 100;
Ny = 1;
CS_flag = 0; % Use complex step, otherwise use centered finite difference
parallel_flag=0;
[G_f,Lambda,W,act_scores] = ActiveSubspace_SVD_Algorithm(f,UB,LB,M,Ny,CS_flag,param_ids,param_base,parallel_flag);
%%

figure(10);clf;
subplot(1,2,1);
bar(diag(Lambda)./max(diag(Lambda)));
xticks(1:length(param_ids));
title('Eigenvalues/Singular Values'); grid on;
set(gca,'FontSize',20);

subplot(1,2,2);
bar(act_scores./max(act_scores))
xticks(1:length(param_ids));
xticklabels(param_names(param_ids));
title('Activity Scores'); grid on;
set(gca,'FontSize',20);

%% See if we can't find a reduced subset capable of predicting the QoI
theta0 = param(param_ids)';
r = 2;
n_train = 500;
act_subspace = @(theta) W(:,1:r)'*theta;


x_train = zeros(r,n_train);
y_train = zeros(1,n_train);


lhs_samp = lhsdesign(n_train,length(param_ids));
par_train = LB(param_ids)+(UB(param_ids)-LB(param_ids)).*lhs_samp;
for j=1:n_train
    par_eval = param_base; par_eval(param_ids) = par_train(j,:);
    y_train(j) = f(par_eval);

    x_train(:,j) = act_subspace(par_train(j,:)');
end
if r==1
    [x_train,sort_ID] = sort(x_train);
    y_train = y_train(sort_ID);
end


X_poly = [ones(1,n_train); x_train; x_train.^2; x_train.^3; x_train.^4];

B = inv(X_poly*X_poly')*X_poly*y_train';
%%
if r==1
    figure(3);clf; hold on;
    % plot(y_train'-X_poly'*B,'ko');
    plot(x_train,y_train,'ko')
    plot(x_train,X_poly'*B,'--r');
    title('OLS Regression')
    grid on; set(gca,'FontSize',20);
end

figure(4);clf; hold on;
plot(y_train'-X_poly'*B,'ko');
title('Residual - Act Subspace Linear Model')
grid on; set(gca,'FontSize',20);

%% Now look at doing MCMC with active subspaces
param_star = param;
UB_full = 1.25.*param_star;
LB_full = 0.75.*param_star;
prior_gauss = @(param,mu,covar) exp(-(param(:)-mu(:))'*inv(covar)*(param(:)-mu(:)))./sqrt(2.*pi.*det(covar));


prior_mu_full = param_star; % Twenty percent above
prior_var_full = 0.05.*prior_mu_full;
prior_covar_full = eye(num_param,num_param).*prior_var_full; % Use a covariance since this is a 3D problem
prior_full = @(param) prior_gauss(param,prior_mu_full,prior_covar_full);


UB_red= max(act_subspace(UB_full),act_subspace(LB_full)); % Twenty percent above
LB_red = min(act_subspace(UB_full),act_subspace(LB_full)); % Twenty percent above

prior_mu_red = act_subspace(param_star); % Twenty percent above
prior_var_red = 0.05.*prior_mu_red;
prior_covar_red = eye(num_param,num_param).*prior_var_red; % Use a covariance since this is a 3D problem

prior_red = @(param) prior_gauss(param,prior_mu_red,prior_covar_red);

%% Chain when covariance is not from sensitivity
actsub_lin = @(q,dummy) sqrt([1; q; q.^2; q.^3; q.^4]'*B);
MC_samp = 10000;
k0=1000;
covar_full = eye(num_param,num_param).*0.01.*abs(param_guess);
covar_red = eye(r,r).*0.01.*abs(act_subspace(param_guess));

%%
tic
[chain_full,s2_full] = ...
AdaptiveMetropolis_s2update(f_mod,data,prior_full,param_guess',k0,MC_samp,covar_full,UB_full,LB_full);
toc
%%
tic
[chain_red,s2_red] =...
AdaptiveMetropolis_s2update(actsub_lin,0.*data,prior_red,act_subspace(param_guess),k0,MC_samp,covar_red,UB_red,LB_red);
toc
%%
figure(1);clf
subplot(2,3,1); hold on;
plot(chain_full(1,:));
yline(param_star(1),'--k','LineWidth',2);
subplot(2,3,2); hold on;
plot(chain_full(2,:));
yline(param_star(2),'--k','LineWidth',2);
subplot(2,3,3); hold on;
plot(chain_full(3,:));
yline(param_star(3),'--k','LineWidth',2);
subplot(2,3,4); hold on;
plot(chain_full(4,:));
yline(param_star(4),'--k','LineWidth',2);
subplot(2,3,5); hold on;
plot(chain_full(5,:));
yline(param_star(5),'--k','LineWidth',2);

figure(20);
plot(chain_red');
yline(act_subspace(param_star),'--k','LineWidth',2);


chain_act = W(:,1:r)*chain_red;
chain_inact = W(:,r+1:end)*normrnd(0,1,num_param-r,MC_samp);
chain_transform = chain_act+chain_inact;
figure(11);clf;
subplot(2,3,1); hold on;
plot(chain_transform (1,:));
yline(param_star(1),'--k','LineWidth',2);
subplot(2,3,2); hold on;
plot(chain_transform (2,:));
yline(param_star(2),'--k','LineWidth',2);
subplot(2,3,3); hold on;
plot(chain_transform (3,:));
yline(param_star(3),'--k','LineWidth',2);
subplot(2,3,4); hold on;
plot(chain_transform (4,:));
yline(param_star(4),'--k','LineWidth',2);
subplot(2,3,5); hold on;
plot(chain_transform (5,:));
yline(param_star(5),'--k','LineWidth',2);

%%
figure(2);clf
subplot(2,3,1); hold on;
ksdensity(chain(1,:));
xline(param_star(1),'--k','LineWidth',2);
subplot(2,3,2); hold on;
ksdensity(chain(2,:));
xline(param_star(2),'--k','LineWidth',2);
subplot(2,3,3); hold on;
ksdensity(chain(3,:));
xline(param_star(3),'--k','LineWidth',2);
subplot(2,3,4); hold on;
ksdensity(chain(4,:));
xline(param_star(4),'--k','LineWidth',2);
subplot(2,3,5); hold on;
ksdensity(chain(5,:));
xline(param_star(5),'--k','LineWidth',2);

%%

%% Plot pairwise and get MAP points
MAP_const = zeros(num_param,1);
MAP_est = zeros(num_param,1);

figure(5); clf;
for i=1:num_param
    for j=i+1:num_param
        subplot(num_param-1,num_param-1,(num_param-1)*(i-1)+j-1); hold on;
        plot(chain(j,:),chain(i,:),'o');
    end
    % Get density values
    [f,x] = ksdensity(chain(i,:));
    % Find the maximum of the PDF
    [who,where] = max(f);
    % Store it
    MAP_const(i) = x(where);

end
%%

disp('MAP - const');
disp(MAP_const');
disp('True');
disp(param_star');

figure(10);clf; hold on;
plot(tspace,f_mod(MAP_const),'LineWidth',2);
plot(tspace,data,'ko');
plot(tspace,true_signal,'--k','LineWidth',2)
legend('MAP','Data','Truth')
%%
function out = act_subspace_f(q,IC,tspace,outflag,data)
param = exp(q);
% param = q;
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
