% Metropolis_Algorithm

function [chain,s2chain] = AdaptiveMetropolis_s2update(f,data,prior_F,theta0,k0,M,covar,UB,LB)

n_par = length(theta0);
chain = zeros(n_par,M);
s2chain = zeros(1,M);
n_y   = length(data);

if isempty(covar)
    covar = eye(n_par).*0.1.*(abs(theta0));
end

ss_func = @(param,data) sum((f(param)-data).^2);
lik_func = @(ss,s2) exp( -ss./(2.*s2) ).* (2.*pi.*s2).^(-n_y/2);

chain(:,1) = theta0;
ss_old = ss_func(theta0,data);
sig_s = ss_old./(n_y-n_par); % Estimate of error variance
lik_old = lik_func(ss_old,sig_s);
prior_old = prior_F(theta0);
R = chol(covar);

%% Error variance hyperparameters
n_s   = 0.01; % For inverse gamma
aval = 0.5.*(n_s+n_y); % This is constant
bval = 0.5.*(n_s.*sig_s+ss_old);
s2chain(1) = 1./gamrnd(aval,1./bval);

%% Initialize covariance estimates
Vk = covar;
theta_bar = theta0(:);
Ip = eye(n_par).*1e-8; %Small jitter on diagonal for stability
sp_cov = 2.38.^2 ./ n_par;
%%
num_acc = 0;

    prob_acc = unifrnd(0,1,M,1);
for i=2:M
    theta_star = chain(:,i-1)+R'*normrnd(0,1,n_par,1);

    if sum(any(theta_star(:)>UB(:)))>0 || sum(any(theta_star(:)<LB(:)))>0
        ss_star  = -inf;
        lik_star = -inf;
    else
        ss_star = ss_func(theta_star,data);
        lik_star = lik_func(ss_star,s2chain(i-1));
        lik_old  = lik_func(ss_old,s2chain(i-1)); % Note: we recompute so that the likelihood is compared for the same s2 value
    end
    prior_star = prior_F(theta_star);

    acc_prob = prior_star.*lik_star./(prior_old.*lik_old);
    if acc_prob>prob_acc(i)
        chain(:,i) = theta_star;
        ss_old = ss_star;
        lik_old = lik_star;
        prior_old = prior_star;
        num_acc = num_acc+1;
    else
        chain(:,i) = chain(:,i-1);
    end

    % Now update sigma squared using an inverse gamma (1/gammarnd)
    bval = 0.5.*(n_s.*sig_s+ss_old);
    s2chain(i) = 1./gamrnd(aval,1./bval);

    % Update values of covariance and mean parameter
    
    theta_bar_old = theta_bar;
    theta_bar = (i./(i+1)).*theta_bar_old+chain(:,i)./(i+1);
    Vk = (i-1)*Vk./i + (sp_cov./i).*(i.*theta_bar_old*theta_bar_old' ...
        - (i+1).*theta_bar*theta_bar' + chain(:,i)*chain(:,i)' + Ip);
    
    if mod(i,k0)==0 % update covariance in sampling
        R = chol(Vk);
    end
end
disp('Acceptance rate')
disp((num_acc./M) .* 100);
end