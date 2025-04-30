%  Morris_Screening_Algorithm
% ORIGINAL AUTHOR: Mitchel J. Colebank
% Script for MATH 728: UQ for Phys and Biol Sys
% Date created: February, 2025
%
% Uses the trajectory algorithm to identify Morris' indices
function [mu,mu_star,sigma] = Morris_Screening_Algorithm(f,UB,LB,M,param_ids,param_base)

%%
ids_fix = 1:length(UB);
ids_fix(param_ids) = [];
num_par = length(param_ids);
UB = UB(param_ids);
LB = LB(param_ids);

%%
p = num_par;   % Number of parameters
l = 60;                % Number of levels
delta = l./(2*(l-1));  % Step Size
upper = UB(:)';
lower = LB(:)';
d = zeros(M,p); % Store the elementary effects
%% Try to use the randomization algorithm
% Note that all parameters are scaled to be in the range 0,1 and then
% rescaled in the model evaluation.
A = zeros(p+1,p);
for i=1:p
    A(i+1:p+1,i) = ones((p-(i-1)),1);
end
X = zeros((p+1).*M,p.*M);
%% Random sampling algorithm with random perturbations
qstar = unifrnd(0,1,M,p);
q_eval = zeros(1,length(param_ids)+length(ids_fix));
q_eval(ids_fix) = param_base(ids_fix);
Jp = ones(p+1,p);
J1 = ones(p+1,1);
P = eye(p,p);
UL_MAT = eye(p).*(upper-lower);
for i=1:M
    qcurr = qstar(i,:);
    pm1 = rand(p,1);
    Dstar = eye(p).*(pm1 > 0.5) - eye(p).*(pm1 <= 0.5);
    [~,where] = sort(rand(p,1));
    Pstar = P(where,:);
    Astar = J1*qcurr + (delta./2).*(( (2.*A - Jp)*Dstar + Jp))*Pstar;
    C = J1*(lower) + Astar*UL_MAT;
    q_eval(param_ids) = C(1,:);
    fpast = f(q_eval);
    for j=1:p
        q_eval(param_ids) = C(j+1,:);
        fstep = f(q_eval);
        % Calculate the elementary effect with the QoI
        par_sign = sign(C(j+1,where(j)) - C(j,where(j)));
        d(i,where(j)) =  par_sign.*(fstep - fpast)./delta; % Elementary effect
        fpast = fstep;
    end
    if mod(i,10)==0
        fprintf("%d percent complete.\n",i/M*100)
    end
end

mu = mean(d,1);
mu_star= mean(abs(d),1);
sigma = sqrt(sum((d-mu_star).^2,1)./(M-1));

end

