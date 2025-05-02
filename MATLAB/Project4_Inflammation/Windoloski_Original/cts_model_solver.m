function sol = cts_model_solver(data)

tspan = data.tspan;
ic = data.ic;

pars = exp(data.pars);

options = odeset('RelTol',1e-8,'AbsTol',1e-8); 
sol = ode15s(@(t,y) cts_model_RHS(t,y,pars),tspan,ic,options);