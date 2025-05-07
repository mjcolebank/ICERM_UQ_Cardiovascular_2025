% ActiveSubspace_SVD_Algorithm.m
% Author: Mitchel J. Colebank
% Script for ICERM 2025 workshop on UQ in Math Bio
% Date created: 5/5/2025
%
% Uses active subspace methods (see Constantine  ISBN 978-1-611973-85-3)
% Outputs:
% G_f: Mxp matrix representing the scaled gradient with respect to the parameters
% Lambda: the singular values of G_f
% W: the rotation matrix of G_f
% act_scores: the activity scores of the model parameters, corresponding to
% the parameter influence and model sensitivity
function [G_f,Lambda,W,act_scores] = ActiveSubspace_SVD_Algorithm(f,UB,LB,M,Ny,CS_flag,param_ids,param_base,parallel_flag)
% %% Construct active subspace representation

ids_fix = 1:length(UB);
ids_fix(param_ids) = [];

par_LHS = zeros(length(param_ids)+length(ids_fix),M);


num_param = length(param_ids);
UB = UB(param_ids); UB = UB(:);
LB = LB(param_ids); LB = LB(:);
X = lhsdesign(M,num_param); % Sample [0,1], then scale back to original input space
h = 1e-6;
par_LHS(param_ids,:) = LB+(UB-LB).*X';
par_LHS(ids_fix,:) = ones(length(ids_fix),M).*param_base(ids_fix)';

par_LHS = par_LHS';
grad_f = zeros(num_param,Ny);
G_f = zeros(num_param,M);

%% No Parallel
if parallel_flag==0
    if CS_flag==0 % Finite difference, no complex step
        I_h = eye(length(param_ids)+length(ids_fix)).*h;
        for i=1:M
            for j=1:num_param
                id = param_ids(j);
                fstepP = f(par_LHS(i,:)+I_h(id,:));
                fstepM = f(par_LHS(i,:)-I_h(id,:));
                grad_f(j,:) =  0.5.*(fstepP-fstepM)./h;
            end
            if mod(i/M*100,10)==0
                fprintf("%d percent complete.\n",i/M*100)
            end
            G_f(:,i) = grad_f;
        end
    else % Use complex step
        I_h = eye(length(param_ids)+length(ids_fix)).*h.*1i;
        for i=1:M
            for j=1:num_param
                id = param_ids(j);
                fstepP = f(par_LHS(i,:)+I_h(id,:));
                grad_f(j,:) =  imag(fstepP)./h;
            end
            if mod(i/M*100,10)==0
                fprintf("%d percent complete.\n",i/M*100)
            end
            G_f(:,i) = grad_f;
        end
    end
%% Parallel processing
else 
    if CS_flag==0 % Finite difference, no complex step
        I_h = eye(length(param_ids)+length(ids_fix)).*h;
        for i=1:M
            parfor j=1:num_param
                id = param_ids(j);
                fstepP = f(par_LHS(i,:)+I_h(id,:));
                fstepM = f(par_LHS(i,:)-I_h(id,:));
                grad_f(j,:) =  0.5.*(fstepP-fstepM)./h;
            end
            G_f(:,i) = grad_f;
        end
    else % Use complex step
        I_h = eye(length(param_ids)+length(ids_fix)).*h.*1i;
        for i=1:M
            parfor j=1:num_param
                id = param_ids(j);
                fstepP = f(par_LHS(i,:)+I_h(id,:));
                grad_f(j,:) =  imag(fstepP)./h;
            end
            if mod(i/M*100,10)==0
                fprintf("%d percent complete.\n",i/M*100)
            end
            G_f(:,i) = grad_f;
        end
    end
end
G_f = G_f./sqrt(M);


[W,Lambda,~] = svd(G_f,'econ');
act_scores = sum(Lambda.^2*W'.^2);

end

%% NOTE: to build active subspace, you can transform the parameter 
% dimensionality to an r<p subspace given by
% act_subspace = W(:,1:r)'*theta
% where theta is the original parameters.



% r = 1;
% n_train = 500;
% act_subspace = @(theta) W(:,1:r)'*theta;
%
%
% x_train = zeros(r,n_train);
% y_train = zeros(1,n_train);
%
%
% lhs_samp = lhsdesign(n_train,num_param);
% par_train = LB+(UB-LB).*lhs_samp';
% for j=1:n_train
%     y_train(j) = f_OLS(par_train(:,j));
%
%     x_train(:,j) = act_subspace(par_train(:,j));
% end
% [x_train,sort_ID] = sort(x_train);
% y_train = y_train(sort_ID);
%
%
% X_poly = [ones(1,n_train); x_train; x_train.^2; x_train.^3; x_train.^4];
%
% B = inv(X_poly*X_poly')*X_poly*y_train';
%
% figure(3);clf; hold on;
% % plot(y_train'-X_poly'*B,'ko');
% plot(x_train,y_train,'ko')
% plot(x_train,X_poly'*B,'--r');
% title('OLS Regression')
% grid on; set(gca,'FontSize',20);
%
% %% Try using the linear model for optimization
% x_val_transform = @(q) [1; act_subspace(q); act_subspace(q).^2; act_subspace(q).^3; act_subspace(q).^4]';
% actsub_OLS = @(q) x_val_transform(q)*B;
% param_guess = param_star.*unifrnd(0.9,1.1,num_param,1);
% [param_opt,J_final] = fminsearch(actsub_OLS,param_guess);
%
% % COmpare to full model
% J_full = @(q) OLS_spring(q,IC,tspace,data);
% [param_opt_full,J_final_full] = fminsearch(J_full,param_guess);
%
% disp([param_opt_full param_opt param_star]);
% disp('error in full vs Act subspace');
% disp([(param_opt_full-param_star)./param_star (param_opt-param_star)./param_star])
% disp(sum([(param_opt_full-param_star)./param_star (param_opt-param_star)./param_star].^2))
%
%
% figure(4); clf; hold on;
% plot(tspace,f_mod(param_opt),'--r');
% plot(tspace,f_mod(param_opt_full),'-.b');
% plot(tspace,true_signal,'k')
% plot(tspace,data,'ko');
%
% %% We can also use Active Subspace emulation to speed up MCMC
% UB_full = 1.25.*param_star;
% LB_full = 0.75.*param_star;
% prior_gauss = @(param,mu,covar) exp(-(param(:)-mu(:))'*inv(covar)*(param(:)-mu(:)))./sqrt(2.*pi.*det(covar));
%
%
% prior_mu_full = param_star; % Twenty percent above
% prior_var_full = 0.05.*prior_mu_full;
% prior_covar_full = eye(num_param,num_param).*prior_var_full; % Use a covariance since this is a 3D problem
% prior_full = @(param) prior_gauss(param,prior_mu_full,prior_covar_full);
%
%
% UB_red= max(act_subspace(UB_full),act_subspace(LB_full)); % Twenty percent above
% LB_red = min(act_subspace(UB_full),act_subspace(LB_full)); % Twenty percent above
%
% prior_mu_red = act_subspace(param_star); % Twenty percent above
% prior_var_red = 0.05.*prior_mu_red;
% prior_covar_red = eye(num_param,num_param).*prior_var_red; % Use a covariance since this is a 3D problem
%
% prior_red = @(param) prior_gauss(param,prior_mu_red,prior_covar_red);
%
% %% Chain when covariance is not from sensitivity
% actsub_lin = @(q,dummy) sqrt([1; q; q.^2; q.^3; q.^4]'*B);
% MC_samp = 10000;
% k0=1000;
% covar_full = eye(num_param,num_param).*0.01.*abs(param_guess);
% covar_red = eye(r,r).*0.01.*abs(act_subspace(param_guess));
%
% %%
% tic
% [chain_full,s2_full] = ...
%     AdaptiveMetropolis_s2update(f_mod,data,prior_full,param_guess',k0,MC_samp,covar_full,UB_full,LB_full);
% toc
% %%
% tic
% [chain_red,s2_red] =...
%     AdaptiveMetropolis_s2update(actsub_lin,0.*data,prior_red,act_subspace(param_guess),k0,MC_samp,covar_red,UB_red,LB_red);
% toc
% %%
% figure(1);clf
% subplot(2,3,1); hold on;
% plot(chain_full(1,:));
% yline(param_star(1),'--k','LineWidth',2);
% subplot(2,3,2); hold on;
% plot(chain_full(2,:));
% yline(param_star(2),'--k','LineWidth',2);
% subplot(2,3,3); hold on;
% plot(chain_full(3,:));
% yline(param_star(3),'--k','LineWidth',2);
% subplot(2,3,4); hold on;
% plot(chain_full(4,:));
% yline(param_star(4),'--k','LineWidth',2);
% subplot(2,3,5); hold on;
% plot(chain_full(5,:));
% yline(param_star(5),'--k','LineWidth',2);
%
% figure(20);
% plot(chain_red');
% yline(act_subspace(param_star),'--k','LineWidth',2);
%
%
% chain_act = W(:,1:r)*chain_red;
% chain_inact = W(:,r+1:end)*normrnd(0,1,num_param-r,MC_samp);
% chain_transform = chain_act+chain_inact;
% figure(11);clf;
% subplot(2,3,1); hold on;
% plot(chain_transform (1,:));
% yline(param_star(1),'--k','LineWidth',2);
% subplot(2,3,2); hold on;
% plot(chain_transform (2,:));
% yline(param_star(2),'--k','LineWidth',2);
% subplot(2,3,3); hold on;
% plot(chain_transform (3,:));
% yline(param_star(3),'--k','LineWidth',2);
% subplot(2,3,4); hold on;
% plot(chain_transform (4,:));
% yline(param_star(4),'--k','LineWidth',2);
% subplot(2,3,5); hold on;
% plot(chain_transform (5,:));
% yline(param_star(5),'--k','LineWidth',2);
%
% %%
% figure(2);clf
% subplot(2,3,1); hold on;
% ksdensity(chain(1,:));
% xline(param_star(1),'--k','LineWidth',2);
% subplot(2,3,2); hold on;
% ksdensity(chain(2,:));
% xline(param_star(2),'--k','LineWidth',2);
% subplot(2,3,3); hold on;
% ksdensity(chain(3,:));
% xline(param_star(3),'--k','LineWidth',2);
% subplot(2,3,4); hold on;
% ksdensity(chain(4,:));
% xline(param_star(4),'--k','LineWidth',2);
% subplot(2,3,5); hold on;
% ksdensity(chain(5,:));
% xline(param_star(5),'--k','LineWidth',2);
%
% %%
%
% %% Plot pairwise and get MAP points
% MAP_const = zeros(num_param,1);
% MAP_est = zeros(num_param,1);
%
% figure(5); clf;
% for i=1:num_param
%     for j=i+1:num_param
%         subplot(num_param-1,num_param-1,(num_param-1)*(i-1)+j-1); hold on;
%         plot(chain(j,:),chain(i,:),'o');
%     end
%     % Get density values
%     [f,x] = ksdensity(chain(i,:));
%     % Find the maximum of the PDF
%     [who,where] = max(f);
%     % Store it
%     MAP_const(i) = x(where);
%
% end
% %%
%
% disp('MAP - const');
% disp(MAP_const');
% disp('True');
% disp(param_star');
%
% figure(10);clf; hold on;
% plot(tspace,f_mod(MAP_const),'LineWidth',2);
% plot(tspace,data,'ko');
% plot(tspace,true_signal,'--k','LineWidth',2)
% legend('MAP','Data','Truth')
% %%
%
%
%
%     function J = OLS_spring(q,IC,tspace,data)
%         y = call_spring_model(q,IC,tspace);
%         J = sum((y(:)-data(:)).^2);
%     end
%
%     function y = call_spring_model(q,IC,tspace)
%         [~,y] = ode45(@spring_model_5par,tspace,IC,[],q);
%         y = y(:,1)';
%     end
%
%     function dy = spring_model_5par(t,y,params)
%
%         % Unpack parameters
%         m      = params(1);
%         c      = params(2);
%         k      = params(3);
%         f0     = params(4);
%         omegaF = params(5);
%
%         % Redefine state variables for convenience
%         x = y(1);
%         v = y(2);
%
%
%         % RHS equations
%         accel_eq = (- c.*v - k.*x)./m + f0.*cos(omegaF.*t)/m;
%         dy = [v;
%             accel_eq];
%
%     end