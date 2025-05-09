% Metropolis_Regression.m
% Author: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: March, 2025
%
clear; clc; close all;
addpath('Project0_UQMethods\')
%% Compute the posterior distribution under the assumption of a uniform
% Initial conditions
z0 = 2; v0 = -2;
IC = [z0;v0];

% Five parameter non-homogeneous spring problem
% The spring model has three parameters
m = 3;        % mass, kg
c = 2;        % damper friction, kg/s
k = 5;        % spring resistance, kg/s^2
F0 = 1;     % forcing term, kg cm / s^2
omegaF = 0.8;%1.8; % frequency of forcing term, rad/s
noise_std = 0.1;

% Put all the parameters together
param_star = [m;c;k;F0;omegaF]; % Get rid of mass term
par_names = {'m','c','k','F0','\omega'};
tend = 4;
%% 
tspace = linspace(0,tend,20);
num_param = length(param_star);
n_xpts = length(tspace);


f_mod = @(q) call_spring_model(q,IC,tspace);
true_signal = f_mod(param_star);
n_data = length(true_signal);
data = true_signal+normrnd(0,noise_std,1,n_data);

f_OLS = @(q) OLS_spring(q,IC,tspace,data);

%% Construct active subspace representation
UB = 1.5.*param_star;
LB = 0.5.*param_star;
M = 2000;
Ny = 1;
CS_flag=1;
param_ids = 1:5;
[G_f,S,W,act_scores] = ActiveSubspace_SVD_Algorithm(f_OLS,UB,LB,M,Ny,CS_flag,param_ids,param_star',0);

%%
figure(10);clf;
subplot(1,2,1);
bar(diag(S)./max(S(:)));
xticks(1:num_param)
title('Eigenvalues'); grid on;
set(gca,'FontSize',20);

subplot(1,2,2);
bar(act_scores./max(act_scores))
xticks(1:num_param)
xticklabels(par_names);
title('Activity Scores'); grid on;
set(gca,'FontSize',20);

%% Now build active subspace
r = 2;
n_train = 500;
% act_subspace = @(theta) W(1:r,:)*theta;
act_subspace = @(theta) W(:,1:r)'*theta;

% Z = (W(1:r,:)*param_star)
% Z = (W(:,1:r)'*param_star)

x_train = zeros(r,n_train);
y_train = zeros(1,n_train);


lhs_samp = lhsdesign(n_train,num_param);
par_train = LB+(UB-LB).*lhs_samp';
for j=1:n_train
    y_train(j) = f_OLS(par_train(:,j));
   
    x_train(:,j) = act_subspace(par_train(:,j));
end
AS_max = max(x_train,[],2);
AS_min = min(x_train,[],2);
% [x_train,sort_ID] = sort(x_train);
% y_train = y_train(sort_ID);


X_poly = [ones(1,n_train); x_train; x_train.^2];%; x_train.^3; x_train.^4];

 B = inv(X_poly*X_poly')*X_poly*y_train';
 
 figure(3);clf; hold on;
 plot(y_train'-X_poly'*B,'ko');
 % plot(x_train,y_train,'ko')
 % plot(x_train,X_poly'*B,'--r');
 title('OLS Regression')
 grid on; set(gca,'FontSize',20);

 %% Try using the linear model for optimization
 x_val_transform = @(q) [1; act_subspace(q); act_subspace(q).^2]';%; act_subspace(q).^3; act_subspace(q).^4]';
 actsub_OLS = @(q) x_val_transform(q)*B;
 param_guess = param_star.*unifrnd(0.7,1.3,num_param,1);
 [param_opt,J_final] = fminsearch(actsub_OLS,param_guess);

 % COmpare to full model
 J_full = @(q) OLS_spring(q,IC,tspace,data);
 [param_opt_full,J_final_full] = fminsearch(J_full,param_guess);

 disp([param_opt_full param_opt param_star]);
 disp('error in full vs Act subspace');
 disp([(param_opt_full-param_star)./param_star (param_opt-param_star)./param_star])
 disp(sum([(param_opt_full-param_star)./param_star (param_opt-param_star)./param_star].^2))


 figure(4); clf; hold on;
 plot(tspace,f_mod(param_opt),'--r');
 plot(tspace,f_mod(param_opt_full),'-.b');
 plot(tspace,true_signal,'k')
 plot(tspace,data,'ko');
 legend('AS','Full','truth','data')
 
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

prior_mu_red = act_subspace(param_star); % Twenty percent above
prior_var_red = 0.01.*abs(prior_mu_red);
prior_covar_red = eye(r,r).*prior_var_red; % Use a covariance since this is a 3D problem

prior_red = @(param) prior_gauss(param,prior_mu_red,prior_covar_red);

%% Chain when covariance is not from sensitivity
% invsubspace = @(q) 
actsub_lin = @(q,dummy) [1; q; q.^2]'*B;
actsub_ODE = @(q) OLS_spring(q,IC,tspace,data);
% actsub_lin = @(q,dummy) sqrt([1; q; q.^2; q.^3; q.^4]'*B);
MC_samp = 10000;
k0=1000;
covar_full = eye(num_param,num_param).*0.01.*abs(param_guess);
covar_red = eye(r,r).*0.001.*abs(act_subspace(param_guess));

%%
tic
[chain_full,s2_full] = ...
AdaptiveMetropolis(f_OLS,data,prior_full,param_guess',k0,MC_samp,covar_full,UB_full,LB_full);
toc
%%
tic
[chain_red,s2_red] =...
AdaptiveMetropolis(actsub_lin,0.*data,prior_red,act_subspace(param_guess),k0,MC_samp,covar_red,UB_red,LB_red);
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
%%
figure(20);
plot(chain_red');
yline(act_subspace(param_star),'--k','LineWidth',2);

%%
chain_act = W(:,1:r)*chain_red;
chain_inact = W(:,r+1:end)*W(:,r+1:end)'*(normrnd(0,1,num_param,MC_samp)+param_guess);
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
ksdensity(chain_full(1,:));
ksdensity(chain_transform(1,:));
xline(param_star(1),'--k','LineWidth',2);

subplot(2,3,2); hold on;
ksdensity(chain_full(2,:));
ksdensity(chain_transform(2,:));
xline(param_star(2),'--k','LineWidth',2);

subplot(2,3,3); hold on;
ksdensity(chain_full(3,:));
ksdensity(chain_transform(3,:));
xline(param_star(3),'--k','LineWidth',2);

subplot(2,3,4); hold on;
ksdensity(chain_full(4,:));
ksdensity(chain_transform(4,:));
xline(param_star(4),'--k','LineWidth',2);

subplot(2,3,5); hold on;
ksdensity(chain_full(5,:));
ksdensity(chain_transform(5,:));
xline(param_star(5),'--k','LineWidth',2);

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



function J = OLS_spring(q,IC,tspace,data)
y = call_spring_model(q,IC,tspace);
J = sum((y(:)-data(:)).^2);
end

function y = call_spring_model(q,IC,tspace)
[~,y] = ode45(@spring_model_5par,tspace,IC,[],q);
y = y(:,1)';
end

function dy = spring_model_5par(t,y,params)

% Unpack parameters
m      = params(1);
c      = params(2);
k      = params(3);
f0     = params(4);
omegaF = params(5);

% Redefine state variables for convenience
x = y(1);
v = y(2);


% RHS equations
accel_eq = (- c.*v - k.*x)./m + f0.*cos(omegaF.*t)/m;
dy = [v;
      accel_eq];

end