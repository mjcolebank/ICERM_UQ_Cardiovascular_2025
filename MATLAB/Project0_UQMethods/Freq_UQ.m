
function [par_CI,Y_CI,Y_PI] = Freq_UQ(f,param_opt,data,CS_flag,param_ids,param_base)
n_xpts = length(data);
num_param = length(param_ids);
ids_fix = 1:length(UB);
ids_fix(param_ids) = [];
param_all = param_base;
param_all(par_ids) = param_opt;
%% After running the optmization, we can construct confidence intervals for the estimates
% and provide output prediction intervals
% First, we need the covariance, sigma^2 * inv(F^T*F)
% We approximate S via finite differences
S = zeros(n_xpts,num_param);
h = 1e-8;
if CS_flag==0
    for i=1:num_param
        param_step = param_opt;
        param_step(i) = param_opt(i)+h;
        S(:,i) = (f(param_step) - f(param_opt))./h;
    end
else
    for i=1:num_param
        param_step = param_opt;
        param_step(i) = param_opt(i)+h.*1i;
        S(:,i) = imag(f(param_step))./h;
    end
end
%%
res = data' - f(param_opt)';
s2 = res'*res./(n_xpts-num_param);
covar = s2.*inv(S'*S);

% Confidence intervals in the parameters
CI_plus  = param_opt'+tinv(0.975,n_xpts-num_param).*sqrt(diag(covar));
CI_minus = param_opt'-tinv(0.975,n_xpts-num_param).*sqrt(diag(covar));

disp('Parameter CI')
disp([CI_minus CI_plus])

%% Now construct response confidence and prediction intervals using a similar formula
% Note that we can extend this to handle "out of sample" test points, e.g.,
% time points in the future or a higher fidelity grid size
% Need model sensitivity at test points
g = S;

Y_pred = f(param_opt);
Y_CI = zeros(n_xpts,2);
Y_PI = zeros(n_xpts,2);
y_stderr_CI = zeros(n_xpts,1);
y_stderr_PI = zeros(n_xpts,1);

for i=1:n_xtest
    y_stderr_CI(i,1) = tinv(0.975,n_xpts-num_param).*sqrt(g(i,:)*covar*g(i,:)');
    Y_CI(i,1) = Y_pred(i) + tinv(0.975,n_xpts-num_param).*sqrt(g(i,:)*covar*g(i,:)');
     Y_CI(i,2) = Y_pred(i) - tinv(0.975,n_xpts-num_param).*sqrt(g(i,:)*covar*g(i,:)');

     y_stderr_PI(i,1) = tinv(0.975,n_xpts-num_param).*sqrt(s2+g(i,:)*covar*g(i,:)');
     Y_PI(i,1) = Y_pred(i) + tinv(0.975,n_xpts-num_param).*sqrt(s2+g(i,:)*covar*g(i,:)');
     Y_PI(i,2) = Y_pred(i) - tinv(0.975,n_xpts-num_param).*sqrt(s2+g(i,:)*covar*g(i,:)');
end

end