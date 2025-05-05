% Solve the solid mechanics of the curved wall patch

function [T,G_f_total,dC_dt, dL_dt,eps_f] = solve_wall(t,Vw,Am,Cm,Am_ref,C,Lsc,pars,AVflag)
% First, determine the wall thickness to cavity radius ratio (dimensionless)
z = 1.5.*(Cm./Am).*Vw;

% Next, calculate the natural myofiber strain
eps_f = 0.5.*log(Am./Am_ref)-(z.^2)./12 - 0.019.*z.^4;

% Now, pass the myocardiac strain to the sarcomere model to obtain relevant
% stress values
[G_f_total,dC_dt, dL_dt] = sarcomere(t,eps_f,Lsc,C,pars,AVflag);


% Finally, calculate the midwall tension.
T = G_f_total.*Vw.*(1.0+z.^2./3 + z.^4./5)./(2.*Am);
end
