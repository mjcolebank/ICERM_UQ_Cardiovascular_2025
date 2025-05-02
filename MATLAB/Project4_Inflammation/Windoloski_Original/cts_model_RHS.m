function xdot = cts_model_RHS(t,y,pars)
                
% ------------------ Model States -------------------

%Inflammation States
E   = y(1);     MR = y(2);      MA   = y(3);      TNF = y(4);
IL6 = y(5);     IL8 = y(6);     IL10 = y(7);

% ------------------- Parameters -------------------

%E, MR, and MA parameters
kE   = pars(1); kMR  = pars(2);  kMA  = pars(3); kMTNF = pars(4); kM = pars(5); 
xME  = pars(6); xM10 = pars(7);  xMTNF = pars(8); 
hME  = pars(9); hM10 = pars(10); hMTNF = pars(11); 
Mmax = pars(12);

%TNF parameters
kTNF   = pars(13); kTNFM = pars(14); 
xTNF10 = pars(15); xTNF6 = pars(16); 
hTNF10 = pars(17); hTNF6 = pars(18); 
wTNF   = pars(19);

%IL6 parameters
k6   = pars(20); k6M = pars(21); k6TNF = pars(22); 
x610 = pars(23); x66 = pars(24); x6TNF = pars(25); 
h610 = pars(26); h66 = pars(27); h6TNF = pars(28); 
w6   = pars(29);

%IL8 parameters
k8   = pars(30); k8M   = pars(31); k8TNF = pars(32); 
x810 = pars(33); x8TNF = pars(34); 
h810 = pars(35); h8TNF = pars(36); 
w8   = pars(37);

%IL10 parameters
k10  = pars(38); k10M = pars(39); k106 = pars(40);
x106 = pars(41); 
h106 = pars(42);
w10  = pars(43);

% ------------- Equations

if (t<=4)
    dE    = 0.5-kE*E;
else
    dE    = -kE*E;
end

dMR   = kMR*MR*(1-MR/Mmax)-upreg(E,xME,hME)*(kM+kMTNF*upreg(TNF,xMTNF,hMTNF))*downreg(IL10,xM10,hM10)*MR;
dMA   = upreg(E,xME,hME)*(kM+kMTNF*upreg(TNF,xMTNF,hMTNF))*downreg(IL10,xM10,hM10)*MR-kMA*MA;
dTNF  = kTNFM*downreg(IL6,xTNF6,hTNF6)*downreg(IL10,xTNF10,hTNF10)*MA-kTNF*(TNF-wTNF);
dIL6  = (k6M+k6TNF*upreg(TNF,x6TNF,h6TNF))*downreg(IL6,x66,h66)*downreg(IL10,x610,h610)*MA-k6*(IL6-w6);
dIL8  = MA*(k8M+k8TNF*upreg(TNF,x8TNF,h8TNF))*downreg(IL10,x810,h810)-k8*(IL8-w8);
dIL10 = (k10M + k106*upreg(IL6,x106,h106))*MA-k10*(IL10 - w10);

xdot = [dE; dMR; dMA; dTNF; dIL6; dIL8; dIL10];

end


% -------------- Hill Functions

function [y] = upreg(state,sat,exponent) %upregulation
  y  = (state)^exponent/(state^exponent + sat^exponent); 
end

function [y] = downreg(state,sat,exponent) %downregulation
  y = sat^exponent/(state^exponent + sat^exponent); 
end


