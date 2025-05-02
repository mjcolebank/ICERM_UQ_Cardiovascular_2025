% Function that defines the Triseg and circulation model
% Returns the differential equations to be solved in the ODE system

function ydot = DE_model_mass(t,y,pars)
% Volumes
Vsa = y(1);
Vsv = y(2);
Vra = y(3);
Vrv = y(4);
Vpa = y(5);
Vpv = y(6);
Vla = y(7);
Vlv = y(8);
Vsw = y(9);
% Sarcomere model states
cla = y(10);
clv = y(11);
cra = y(12);
crv = y(13);
csw = y(14);

Lla = y(15);
Llv = y(16);
Lra = y(17);
Lrv = y(18);
Lsw = y(19);

% y free wall
ysw = y(20);





%% Pericardium
% Calculate the pericardial pressure
V0_heart = pars.Peri(1);
s        = pars.Peri(2);
V_heart = Vla+Vlv+Vra+Vrv;%+Vsw;


p_peri = exp(s.*(V_heart./V0_heart - 1)); % From Jezek



%% Load in parameters
% Parameters for heart chambers
Vpars  = pars.V;
Vwall  = Vpars(1:5);
Am_ref = Vpars(6:10);

V_LVwall = Vwall(2);
V_RVwall = Vwall(4);
V_SWwall = Vwall(5);

% Parameters for the CV models
CVpars = pars.CV;

%%
% Initialize with previous estimates
ym    = ysw; % estimate of the midwall junction radius
Vm_SW = Vsw; %Septal wall volume estimate


%% First, solve equations 9-11 to obtain xm, Am, and Cm
% Note, in original Lumen's code we assume that ysw = ylv = yrv; We do a
% second iteration and update eveything else

Vm_LV = - Vlv - 0.5.*(V_LVwall+V_SWwall)+Vm_SW;
Vm_RV =   Vrv + 0.5.*(V_RVwall+V_SWwall)+Vm_SW;

Vm         = [Vm_LV,Vm_RV,Vm_SW];
c          = [clv,crv,csw];
L          = [Llv,Lrv,Lsw];
Vwall_ven  = Vwall([2 4 5])';
Am_ref_ven = Am_ref([2 4 5])';

%% Now that we have the values of Am and Cm, we want to calculate the
%  mechanics of the curved cardiac walls
[xm,Am,Cm] = get_wallsegment(Vm,ym);
% LV
[Tm_LV,sig_f_LV,c_LV,L_LV,eps_f_LV] = solve_wall(t,Vwall_ven(1),Am(1),Cm(1),Am_ref_ven(1),c(1),L(1),pars,-1);
% RV
[Tm_RV,sig_f_RV,c_RV,L_RV,eps_f_RV] = solve_wall(t,Vwall_ven(2),Am(2),Cm(2),Am_ref_ven(2),c(2),L(2),pars,0);
% SW
[Tm_S,sig_f_S,c_S,L_S,eps_f_S]      = solve_wall(t,Vwall_ven(3),Am(3),Cm(3),Am_ref_ven(3),c(3),L(3),pars,0);

Tm    = [Tm_LV Tm_RV Tm_S];
G_f = [sig_f_LV sig_f_RV sig_f_S];
L     = [L_LV L_RV L_S];
eps_f = [eps_f_LV eps_f_RV eps_f_S];
c     = [c_LV c_RV c_S];
[Tx,Ty] = get_subtensions(Tm,xm,ym);



%%
p_trans_LV = 2.*(Tx(1)./ym);
p_trans_RV = 2.*(Tx(2)./ym);

plv = -p_trans_LV;
prv =  p_trans_RV;\


%% Now calculate variables in the atria
% Calculate midwall volume (see Lumen's code)
Vm_la = 0.5*Vwall(1)+Vla;
Vm_ra = 0.5*Vwall(3)+Vra;

% Calulate atrial curvature
Cm_la = (4.*pi./3./Vm_la).^(1/3);
Cm_ra = (4.*pi./3./Vm_ra).^(1/3);

% Now determine the area of the midwall
Am_la = 4.*pi./(Cm_la.^2);
Am_ra = 4.*pi./(Cm_ra.^2);

% Now get the tensions 
atr_offset = pars.SarcA(end);
[Tm_la,G_f_la,cla,Lla,eps_f_la] = solve_wall(t+atr_offset,Vwall(1),Am_la,Cm_la,Am_ref(1),cla,Lla,pars,1);
[Tm_ra,G_f_ra,cra,Lra,eps_f_ra] = solve_wall(t+atr_offset,Vwall(3),Am_ra,Cm_ra,Am_ref(3),cra,Lra,pars,1);

% Finally, determine the transmural pressure (different than ventricles,
% see Lumens 2009 code)
pla = abs(2.*Tm_la*Cm_la);
pra = abs(2.*Tm_ra*Cm_ra);

%% Account for pericardium
pla = pla+p_peri;
plv = plv+p_peri;
pra = pra+p_peri;
prv = prv+p_peri;


%% Now, consider the rest of the circulatory model
Ra_val = CVpars(1);
Rm_val = CVpars(2);
Rp_val = CVpars(3);
Rt_val = CVpars(4);
Rvc    = CVpars(5);
Rpv    = CVpars(6);
Rs     = CVpars(7);
Rp     = CVpars(8);

Csa    = CVpars(9);
Csv    = CVpars(10);
Cpa    = CVpars(11);
Cpv    = CVpars(12);

psa = Vsa/Csa; % systemic arteries 
psv = Vsv/Csv; % systemic veins 
ppa = Vpa/Cpa; % pressure pulmonary artery
ppv = Vpv/Cpv; % pressure pulmonary vein

% right heart valves
qt_val = max((pra-prv)./Rt_val,0); % tricuspid valve opens 
qp_val = max((prv-ppa)./Rp_val,0); % pulmonic valve opens

% left heart valves
qm_val = max((pla-plv)./Rm_val,0); % mitral valve opens
qa_val = max((plv-psa)./Ra_val,0); % aortic valve opens

% venous flow
qvc = max((psv-pra)./Rvc,0);
qpv = (ppv-pla)./Rpv;

% Non-valve systemic/pulmonary flows
qs   = (psa - psv)/Rs;
qp   = (ppa - ppv)/Rp;



%% dVdt
dVsa = qa_val - qs;
dVsv = qs     - qvc;
dVra = qvc    - qt_val;
dVrv = qt_val - qp_val;
dVpa = qp_val - qp;
dVpv = qp     - qpv;
dVla = qpv    - qm_val;
dVlv = qm_val - qa_val;

dCla = cla;
dClv = c(1);
dCra = cra;
dCrv = c(2);
dCsw = c(3);

dLla = Lla;
dLlv = L(1);
dLra = Lra;
dLrv = L(2);
dLsw = L(3);



%% Differential Equations
ydot        = zeros(20,1);
ydot(1:8)   = [dVsa dVsv dVra dVrv dVpa dVpv dVla dVlv];
ydot(9)     = sum(Tx); %DAE
ydot(10:14) = [dCla dClv dCra dCrv dCsw];
ydot(15:19) = [dLla dLlv dLra dLrv dLsw];
ydot(20)    = sum(Ty); %DAE

end