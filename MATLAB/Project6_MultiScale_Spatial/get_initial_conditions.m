% Compute the Jacobian for the tensions and determine the values of ym and
% Vm_septum necessary for solving the nonlinear root finding problem
% Author: MJ Colebank
% Last edited: 3/10/2022
%
% THIS CODE COMES WITH NO GUARANTEES
function [ym,Vm,Tx0,Ty0] = get_initial_conditions(t,xm,Vm,ym,y,pars,Vv)
c = y([11 13 14])';
L = y([16 18 19])';

Vpars  = pars.V;
Vwall  = Vpars([2 4 5])';
Am_ref = Vpars([7 9 10])';

ym_init = ym;
Vm_sw_init = Vm(3); 

Vm(1) = - Vv(1) - 0.5.*(Vwall(1)+Vwall(3))+Vm(3);
Vm(2) =   Vv(2) + 0.5.*(Vwall(2)+Vwall(3))+Vm(3);

%% Now solve the root finding problem %20 works!!!!! so does 10
options = optimoptions('fsolve','Display','off','MaxIterations',200,'FiniteDifferenceType','central','OptimalityTolerance',1e-16);
F = @(X) tension_roots(X,t,c,L,pars,Vv);
X0 = [ym Vm];

[Xopt,Fval,exitflag,output,J] = fsolve(F,X0,options);
ym = Xopt(1);
Vm = Xopt(2:4);


[xm,Am,Cm] = get_wallsegment(Vm,ym);
[T,~,~,~] = solve_wall(t,Vwall,Am,Cm,Am_ref,c,L,pars,0);
[Tx0,Ty0] = get_subtensions(T,xm,ym);
Tx0 = sum(Tx0);
Ty0 = sum(Ty0);

end

function cost = tension_roots(X,t,c,L,pars,Vv)
Vpars     = pars.V;
Vwall  = Vpars([2 4 5])';
Am_ref = Vpars([7 9 10])';

Vlv = Vv(1);
Vrv = Vv(2);

ym   = X(1);
VmLV = X(2);
VmRV = X(3);
VmSW = X(4);

Vm = X(2:4);

[xm,Am,Cm] = get_wallsegment(Vm,ym);
[T,~,~,~] = solve_wall(t,Vwall,Am,Cm,Am_ref,c,L,pars,0);
[Tx0,Ty0] = get_subtensions(T,xm,ym);
Tx0 = sum(Tx0);
Ty0 = sum(Ty0);

cost = zeros(1,4);
cost(1) = Tx0;
cost(2) = -VmLV - Vlv - 0.5.*(Vwall(1)+Vwall(3))+VmSW;
cost(3) = -VmRV + Vrv + 0.5.*(Vwall(2)+Vwall(3))+VmSW;
cost(4) = Ty0;

end
