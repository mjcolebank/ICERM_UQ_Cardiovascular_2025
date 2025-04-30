%
% Generates the profile likelihood using the nonidentifiable, three
% parameter models
%
% ORIGINALLY WRITTEN BY MJ COLEBANK

clear all
close all

%
t = 0:.1:1;
f = @(theta) theta(1).*theta(2).*t + theta(3);

theta_true = [2 2 3];

output_true = f(theta_true);

theta_bounds = [0 5];
theta_prof = theta_bounds(1):0.1:theta_bounds(2);
n_prof = length(theta_prof);
num_par = 3;
par_ids = 1:num_par;
par_opt = zeros(num_par,num_par,n_prof);
J_PL    = zeros(num_par,n_prof);
for i=1:num_par
    id_infer = par_ids;
    id_fix = par_ids(i);
    id_infer(i) = [];
    par_infer = theta_true(id_infer);
    for ii=1:n_prof
        par_fix = theta_prof(ii);
        [Xopt,FVAL] = fminsearch(@call_PL,par_infer,[],id_infer,par_fix,id_fix,f,output_true);
        par_opt(i,id_infer,ii) = Xopt;
        par_opt(i,id_fix,ii)   = par_fix;
        J_PL(i,ii)             = FVAL;
        par_infer = Xopt;

    end
end


function J = call_PL(par_infer,id_infer,par_fix,id_fix,f,data)
theta = zeros(length(par_infer)+1,1);
theta(id_infer) = par_infer;
theta(id_fix) = par_fix;

J = sum( (f(theta)-data).^2);
end
     
