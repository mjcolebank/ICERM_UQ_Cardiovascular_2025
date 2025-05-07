function xdot = modelBasic(t,y,pars,ts,T)
global Rav Rmv 

%states (volume)
Vao    = y(1);
Vsa    = y(2);
Vsv    = y(3);
Vvc    = y(4);
Vlv    = y(5);

%parameters
Ra  = pars(1); 
Rs  = pars(2); 
Rv	= pars(3);

Eao = pars(4);
Esa = pars(5);
Esv = pars(6);
Evc = pars(7);

Tsf = pars(8);
Trf = pars(9);
Em  = pars(10); 
EMm = pars(11);

%pressures
pao  = Vao*Eao;
psa  = Vsa*Esa;
psv  = Vsv*Esv;
pvc  = Vvc*Evc;
Elv  = ElastanceBasic(t-ts,T,Tsf,Trf,Em,EMm);
plv  = Elv*Vlv; 

%Ohm's law derivd flows - with diodes to represent valves
if pvc > plv
    qvc = (pvc - plv)./Rmv;
else
    qvc = 0;
end
if plv > pao
    qao = (plv - pao)./Rav;
else
    qao = 0;
end
% Flows from Ohms Law
qs   = (psa  - psv)/Rs;
qv   = (psv - pvc)/Rv;
qa   = (pao - psa)/Ra;

% Differential Equations
dVao = qao - qa;
dVsa = qa - qs;
dVsv = qs - qv;
dVvc = qv - qvc;
dVlv = qvc - qao;

xdot = [dVao;...
        dVsa;...
        dVsv;...
        dVvc;...
        dVlv;];            
       
   
        

     
  