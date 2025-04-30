% Wrapper function to get solutions from CV model

function [output,tplot] = call_CV_model(IC,param,tspace)
options = odeset('RelTol',1e-8,'AbsTol',1e-8); % Necessary for the ODE solver

T = param(end);
% Now solve our system of ODEs
y = ode45(@cardiovascular_model,[tspace(1), tspace(end)],IC,options,param);


% We solved the model for 30 heartbeats; we want to extract the last two
% for plotting
tplot = linspace(tspace(end)-T,tspace(end),50);
yout = deval(y,tplot);
tplot = tplot-tplot(1);


% Extract the solutions to the two differential equations
Vlv = yout(1,:);
Vao = yout(2,:);
% param = [P_LA,P_SysCap,Rmv,Rav,Rart,Emax,Emin,T_peak,T_relax,T,Vlv_d,Cao];
P_LA     = param(1);
P_SysCap = param(2);
Rmv      = param(3);
Rav      = param(4);
Rart     = param(5);
Cao      = param(6);
% Now, recompute the pressures and flows
plv = LinearElastance(Vlv,tplot,param(7:12));
pao = Vao./Cao;

% Use the 'max' operator to keep positive flows for valves
qmv = max((P_LA-plv)./Rmv,0);
qav = max((plv-pao)./Rav,0);
qart = (pao-P_SysCap)./Rart;

output.Vlv = Vlv;
output.Vao = Vao;
output.plv = plv;
output.pao = pao;
output.qmv = qmv;
output.qav = qav;
output.qart = qart;

end