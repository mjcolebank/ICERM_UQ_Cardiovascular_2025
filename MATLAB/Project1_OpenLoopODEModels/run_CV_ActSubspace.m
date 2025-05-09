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
load('CV_test_data.mat','data','param_true')
% Starting and ending time
tstart = 0;
tend   = 30.*T; % 30 cycles
dt = 1e-3;
tspace = tstart:dt:tend;
% param_base = log(param);
param_base = param;
%% Run active subspaces depending on the quantity of interest
output_flag = 3;
f = @(q) act_subspace_f(q,IC,tspace,output_flag,data);
param_ids = [1:11];%1:11; % Exclude the cardiac cycle length;
%%
UB = param_base.*1.5;
LB = param_base.*0.5;
% LB(4) = param(4).*0.95;
% UB(4) = param(4).*1.05;
M = 200;
Ny = 1;
CS_flag = 1; % Use complex step, otherwise use centered finite difference
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
param_star = param_true(param_ids)';
r = 1;
n_train = 100;
act_subspace = @(theta) W(:,1:r)'*theta;
inact_subspace = @(theta) W(:,r+1:end)'*theta;


x_train = zeros(r,n_train);
y_train = zeros(1,n_train);


lhs_samp = lhsdesign(n_train,length(param_ids));
% lhs_samp = unifrnd(0,1,n_train,length(param_ids));
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

AS_max = max(x_train,[],2);
AS_min = min(x_train,[],2);
X_poly = [ones(1,n_train); x_train; x_train.^2 ; x_train.^3 ; x_train.^4];

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

rel_res = (y_train'-X_poly'*B);%./y_train';
figure(4);clf; 
plot(rel_res,'ko');
title('Residual - Act Subspace Linear Model')
grid on; set(gca,'FontSize',20);

figure(99);
for i=1:length(param_ids)
subplot(3,4,i); ksdensity(abs(G_f(i,:)));
xlabel(param_names{i});
grid on; set(gca,'FontSize',20);
end
subplot(3,4,2); title('Density of Gradient Values')

 %% Full model optimization
  param_guess = param_star.*unifrnd(0.8,1.2,num_param,1);
 J_full = @(q) act_subspace_sse(q,IC,tspace,output_flag,data,param,param_ids);
 [param_opt_full,J_final_full] = fminsearch(J_full,param_guess);
 %% Try using the linear model for optimization
   param_guess = param_star.*unifrnd(0.8,1.2,num_param,1);

num_param = length(param_ids);
 x_val_transform = @(q) [1; act_subspace(q); act_subspace(q).^2 ; act_subspace(q).^3 ; act_subspace(q).^4]';
 actsub_OLS = @(q) abs(x_val_transform(q)*B);
   [param_opt_red,J_final] = fmincon(actsub_OLS,param_guess,[],[],[],[],LB(param_ids),UB(param_ids));

 % actsub_lin = @(q,dummy) abs([1; q; q.^2]'*B);
 % param_red_guess = act_subspace(param_guess)
 % [param_opt_red,J_final] = fmincon(actsub_lin,param_red_guess,[],[],[],[],AS_min,AS_max);
 % param_opt_red = param_opt_red + W(:,r+1:end)*inact_subspace(param_opt_red)
 % param_opt_red = W(:,1:r)*param_opt_red + W(:,r+1:end)*W(:,r+1:end)'*param_guess
 %%
 disp([param_opt_full param_opt_red param_star]);
 disp('error in full vs Act subspace');
 disp([(param_opt_full-param_star)./param_star (param_opt_red-param_star)./param_star])
 disp(sum([(param_opt_full-param_star)./param_star (param_opt_red-param_star)./param_star].^2))

%%
param_eval_full = param; param_eval_full(param_ids) = param_opt_full;
[f_out_full, ~] = call_CV_model(IC,param_eval_full,tspace);

param_eval_red = param; param_eval_red(param_ids) = param_opt_red;
[f_out_red, ~] = call_CV_model(IC,param_eval_red,tspace);

 figure(4); clf; hold on;
 plot(f_out_full.pao,'--r');
 plot(f_out_red.pao,'-.b');
 plot(data.pao,'ko');
 legend('Full','AS','truth','data')
%% We can also use Active Subspace emulation to speed up MCMC
UB_full = 1.25.*param_star;
LB_full = 0.75.*param_star;
prior_gauss = @(param,mu,covar) exp(-(param(:)-mu(:))'*inv(covar)*(param(:)-mu(:)))./sqrt(2.*pi.*det(covar));


prior_mu_full = param_star; % Twenty percent above
prior_var_full = 0.01.*prior_mu_full;
prior_covar_full = eye(num_param,num_param).*prior_var_full; % Use a covariance since this is a 3D problem
prior_full = @(param) prior_gauss(param,prior_mu_full,prior_covar_full);


UB_red=  AS_max;%max(act_subspace(UB_full),act_subspace(LB_full)); % Twenty percent above
LB_red = AS_min;%(act_subspace(UB_full),act_subspace(LB_full)); % Twenty percent above

prior_mu_red = act_subspace(param_opt_red); % Twenty percent above
prior_var_red = 0.01.*abs(prior_mu_red);
prior_covar_red = eye(r,r).*prior_var_red; % Use a covariance since this is a 3D problem

prior_red = @(param) prior_gauss(param,prior_mu_red,prior_covar_red);

%% Chain when covariance is not from sensitivity
% invsubspace = @(q) 
actsub_lin = @(q,dummy) [1; q; q.^2]'*B;
% actsub_ODE = @(q) OLS_spring(q,IC,tspace,data);
% actsub_lin = @(q,dummy) sqrt([1; q; q.^2; q.^3; q.^4]'*B);
MC_samp = 10000;
k0=1000;
covar_full = eye(num_param,num_param).*0.01.*abs(param_guess);
covar_red = eye(r,r).*0.01.*abs(act_subspace(param_guess));

%%
f_mcmc = @(q) act_subspace_sse(q,IC,tspace,output_flag,data,param,param_ids);
tic
[chain_full,s2_full] = ...
AdaptiveMetropolis(f_mcmc,data,prior_full,param_guess',k0,MC_samp,covar_full,UB_full,LB_full);
toc
%%
tic
[chain_red,s2_red] =...
AdaptiveMetropolis(actsub_lin,data.pao,prior_red,act_subspace(param_guess),k0,MC_samp,covar_red,UB_red,LB_red);
toc
%%
figure(1);clf
for i=1:11
subplot(3,4,i); hold on;
plot(chain_full(i,:));
yline(param_star(i),'--k','LineWidth',2);
end
%%
figure(20);
plot(chain_red');
yline(act_subspace(param_star),'--k','LineWidth',2);

%%

prior_mu_inact = inact_subspace(param_opt_red);
prior_var_inact = 0.01.*abs(prior_mu_inact);
inact_chain = zeros(num_param-r,MC_samp);
for i=1:length(prior_mu_inact)
    inact_chain(i,:) = normrnd(prior_mu_inact(i),sqrt(prior_var_inact(i)),1,MC_samp);
end
%%
chain_act = W(:,1:r)*chain_red;
chain_inact = W(:,r+1:end)*inact_chain;
% chain_inact = sqrt(prior_var_red).*(W(:,r+1:end)*W(:,r+1:end)'*normrnd(0,1,num_param,MC_samp))+param_guess;

chain_transform = chain_act+chain_inact;
figure(11);clf
for i=1:11
subplot(3,4,i); hold on;
plot(chain_transform(i,:));
yline(param_star(i),'--k','LineWidth',2);
end

%%
figure(2);clf
for i=1:11
subplot(3,4,i); hold on;
ksdensity(chain_full(i,:));
ksdensity(chain_transform(i,:));
xline(param_star(i),'--k','LineWidth',2);
xlabel(param_names{param_ids(i)})
end


%%

%% Plot pairwise and get MAP points
MAP_const = zeros(num_param,1);
MAP_est = zeros(num_param,1);

figure(50); clf;
for i=1:num_param
    for j=i+1:num_param
        subplot(num_param-1,num_param-1,(num_param-1)*(i-1)+j-1); hold on;
        plot(chain_full(j,:),chain_full(i,:),'o');
    end
    % Get density values
    [f,x] = ksdensity(chain_full(i,:));
    % Find the maximum of the PDF
    [who,where] = max(f);
    % Store it
    MAP_const(i) = x(where);

end
%%
figure(60); clf;
for i=1:num_param
    for j=i+1:num_param
        subplot(num_param-1,num_param-1,(num_param-1)*(i-1)+j-1); hold on;
        plot(chain_transform(j,:),chain_transform(i,:),'o');
    end
    % Get density values
    [f,x] = ksdensity(chain_transform(i,:));
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

figure(100);clf; hold on;
plot(tspace,f_mod(MAP_const),'LineWidth',2);
plot(tspace,data,'ko');
plot(tspace,true_signal,'--k','LineWidth',2)
legend('MAP','Data','Truth')
%%
function out = act_subspace_f(q,IC,tspace,outflag,data)
% param = exp(q);
param = q;
[output,~] = call_CV_model(IC,param,tspace);

if outflag==1 % systolic LV pressure
    out = max(output.plv);
elseif outflag==2 %systolic Ao pressure
    out = max(output.pao);
elseif outflag==3 %Ao pressure
    res1 = (data.pao - output.pao);
    out = sum(res1.^2);
elseif outflag==4 %residual function (assume pressure in the Aorta, aortic and mitral flow
    res1 = (data.pao - output.pao)./max(data.pao);
    res2 = (data.qav - output.qav)./max(data.qav);
    res3 = (data.qmv - output.qmv)./max(data.qmv);
    out = sum([res1 res2 res3].^2);
end

end

%%
function out = act_subspace_sse(q,IC,tspace,outflag,data,param_all,param_ids)
% param = exp(q);
param = param_all;
param(param_ids) = q;
[output,~] = call_CV_model(IC,param,tspace);

if outflag==1 % systolic LV pressure
    out = max(output.plv);
elseif outflag==2 %systolic Ao pressure
    out = max(output.pao);
elseif outflag==3 %Ao pressure
    res1 = (data.pao - output.pao)./max(data.pao);
    out = sum(res1.^2);
elseif outflag==4 %residual function (assume pressure in the Aorta, aortic and mitral flow
    res1 = (data.pao - output.pao)./max(data.pao);
    res2 = (data.qav - output.qav)./max(data.qav);
    res3 = (data.qmv - output.qmv)./max(data.qmv);
    out = sum([res1 res2 res3].^2);
end

end
