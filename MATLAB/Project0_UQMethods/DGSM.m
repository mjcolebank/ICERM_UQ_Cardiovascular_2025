% DGSM.m
% Author: Mitchel J. Colebank
% Script for ICERM 2025 workshop on UQ in Math Bio
% Date created: 5/5/2025
%
% Computes derivative based global sensitivity metrics
% Outputs:
% mu: the average derivative of the QoI w.r.t. the parameters
% mu_star: the average absolute value of the derivative QoI w.r.t. parameters
% v: the average squared deriative of QoI w.r.t. the parameters
%%
function [mu,mu_star,v] = DGSM(f,UB,LB,M,param_ids,param_base,parallel_flag,CS_flag)
%%
ids_fix = 1:length(UB);
ids_fix(param_ids) = [];
p = length(param_ids);
p_fix = length(ids_fix);
UB = UB(param_ids);
LB = LB(param_ids);
delta = 1e-6;  % Step Size
d = zeros(M,p); % Store the derivatives
%%
mean_value = 0;
mean_sq_value = 0;
% Generate samples for DGSM
qcurr = zeros(M,p+p_fix);
for j=1:p
    qcurr(:,param_ids(j)) = unifrnd(LB(j)+delta,UB(j)-delta,1,M);
end
qcurr(:,ids_fix) = param_base(ids_fix);
%%
if parallel_flag==1 && CS_flag==0 % Parallel with fwd finite difference
    parfor i=1:M
        fbase = f(qcurr(i,:));
        mean_value = mean_value+fbase;
        mean_sq_value = mean_sq_value+fbase.^2;
        for j=1:p
            par_plus = qcurr(i,:);
            par_plus(param_ids(j))=par_plus(param_ids(j))+delta;
            fstep = f(par_plus);
            % Calculate the DGSM
            d(i,j) =  (fstep - fbase)./delta;
            mean_value = mean_value+fstep;
            mean_sq_value = mean_sq_value+fstep.^2;
        end
    end
    %%
elseif parallel_flag==1 && CS_flag==1 % Parallel with complex step
    deltai = delta.*1i;
    parfor i=1:M
        for j=1:p
            par_plus = qcurr(i,:);
            par_plus(param_ids(j))=par_plus(param_ids(j))+deltai;
            fstep = f(par_plus);
            % Calculate the DGSM
            d(i,j) =  imag(fstep)./delta;
            mean_value = mean_value+fstep;
            mean_sq_value = mean_sq_value+fstep.^2;
        end
    end
    %%
elseif parallel_flag==0 && CS_flag==0 % No parallel with Fwd finite difference
    for i=1:M
        fbase = f(qcurr(i,:));
        mean_value = mean_value+fbase;
        mean_sq_value = mean_sq_value+fbase.^2;
        for j=1:p
            par_plus = qcurr(i,:);
            par_plus(param_ids(j))=par_plus(param_ids(j))+delta;
            fstep = f(par_plus);
            % Calculate the DGSM
            d(i,j) =  (fstep - fbase)./delta;
            mean_value = mean_value+fstep;
            mean_sq_value = mean_sq_value+fstep.^2;
        end
    end
    %%
elseif parallel_flag==0 && CS_flag==1 % Parallel with complex step
    deltai = delta.*1i;
    for i=1:M
        for j=1:p
            par_plus = qcurr(i,:);
            par_plus(param_ids(j))=par_plus(param_ids(j))+deltai;
            fstep = f(par_plus);
            % Calculate the DGSM
            d(i,j) =  imag(fstep)./delta;
            mean_value = mean_value+fstep;
            mean_sq_value = mean_sq_value+fstep.^2;
        end
    end
end
%% Calculate DGSM metrics

mu = sum(d,1)./M;
mu_star= sum(abs(d),1)./M;
v = sum(d.^2,1)./M;

%% If you are interested in comparing to Sobol' indices, the DGSM
% is analytically an upper bound on the total Sobol index
% uncomment if wanted
% Assume parameters are uniform, so poincare constant is (b-a).^2./pi.^2
% E[X^2] - E[X]^2
% Calculate the variance of all the outputs
% mean_value = mean_value./(M*(p+1));
% mean_sq_value = mean_sq_value./(M*(p+1));
% var_out = mean_sq_value - mean_value.^2;
% poincare_const = (UB-LB).^2 ./ (pi.^2);
% figure;clf;
% bar(poincare_const.*v./var_out);
% xticklabels(parameter_names);
% set(gca,'FontSize',20);
% grid on;
end

