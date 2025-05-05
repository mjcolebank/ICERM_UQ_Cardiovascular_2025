% Function that takes in the the sarcomere fiber strain and computes the
% total stress as a function of both active and pass stresses.

function [G_f_total,dC_dt,dLsc_dt] = sarcomere(t,eps_f,Lsc,C,pars,AV)
active_scale = 1.0;
if AV==1
    Spars = pars.SarcA;
elseif AV==-1
    Spars = pars.SarcV;
    if length(pars.MI)>1
        active_scale = induce_MI(pars,t);%max(Spars(end-1)*0.1,Spars(end)); %take max of 30% active or the passive
    else
        active_scale = pars.MI;
    end
else
    Spars = pars.SarcV;
end
T     = pars.T;
tmod = mod(t,T);


Ls_ref  = Spars(1);
Ls_iso  = Spars(2);
Vmax    = Spars(3);
Lsc0    = Spars(4);
C_rest  = Spars(5);

% Timing parameters are treated as parameters here: published code uses
% timing relative to some activation time
tauR  = Spars(6);%./active_scale;
tauD  = Spars(7);
tauSC = Spars(8);

% Active and passive fiber stress values
sig_act    = Spars(9);
Ls_ref_pas = Spars(10);
sig_pas    = Spars(11);
k_pas      = Spars(12);


% First, calculate the sarcomere length as a function of Sarcomere strain
Ls = Ls_ref.*exp(eps_f); %eq. B1

% Next, define the differential equation for the contractile element
% length, Lsc
dLsc_dt = ((Ls-Lsc)./Ls_iso - 1).*Vmax; %eq. B2

%% Active force
% Now, define the terms needed for the differential equation for sarcomere
% mechanical activiation, C

CL = tanh(4.0.*(Lsc-Lsc0).^2); %eq. B4

xF = min(8,max(0,tmod/tauR)); % eq. B5b
F_rise = active_scale.*(0.02.*xF.^3 .*(8-xF).^2 .* exp(-xF)); %eq. B5

T_Lsc = tauSC.*(0.29+0.3.*(Lsc)); %eq. B6

% NOTE: This function has a discontinuity at t=T. Lumens actually uses an
% approximation (not included in the text;
exp_func = 1.0+exp((T_Lsc - tmod)./tauD); %in eq. B3
dC_dt = CL.*F_rise./tauR + ((C_rest-C)./exp_func)./tauD; %eq. B3
%

% Calculate the active stress
G_act  = sig_act.*C.*(Lsc-Lsc0).*(Ls-Lsc)./Ls_iso;

%% Passive force (titin and ECM)

%Stretch in the passive components of the wall
Ls_pas = Ls./Ls_ref_pas;

% Nonlinear passive stiffening
G_pas = max(sig_pas.*(Ls_pas.^k_pas - 1.0),0);


G_f_total = G_pas+G_act;

end