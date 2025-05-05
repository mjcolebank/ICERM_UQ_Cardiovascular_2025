function ic = cts_model_ic

global dataC

pars = exp(dataC.pars);

wTNF = pars(19);
w6   = pars(29);
w8   = pars(37);
w10  = pars(43);

%Inflammation initial conditions
TNF0  = wTNF;   
IL60  = w6; 
IL80  = w8;
IL100 = w10;           
MA0   = 0.0;
MR0   = 28200;   
E0    = 0;

% Parameter vector of initial conditions
ic = [E0   MR0   MA0   TNF0   IL60  IL80  IL100];   
end