% Function that provides the nominal parameter values for the
% Triseg/circulation model.
% Author: MJ Colebank
% Last edited: 3/10/2022
%
% THIS CODE COMES WITH NO GUARANTEES
% Abbreviations: LA -> left atrium, LV -> left ventricle, RA -> right
% atrium, RV -> right ventricle, S -> septum, Am -> midwall area, Cm ->
% midwall curvature

function [IC,pars] = get_pars_GR %make variable input in future versions
% Set the nominal parameter values
T = 1.0;

%% Get the average heart rate

SV = 90;
HR = 60;
CO_microl_min = SV*HR;
CO_ml = CO_microl_min/60;




Vtot = 5000;
%% Define conversion factors
p_conv = 0.133322;% mmHg -> Kpa

%Cardiac output: assume 5 L/min
%% First, CV model parameters
% Most are from Boron and Boulapep, Tables 19-3 (vascular) and 22-3 (heart)

p_lv_sys  = 120*p_conv;
p_lv_dias = 5*p_conv;

p_rv_sys  = 20*p_conv;
p_rv_dias = 2*p_conv;

p_sa_sys  = 0.99*p_lv_sys;
p_sa_dias = 0.66*p_lv_sys;%min(p_SA).*p_conv;%
p_sa_mean = (p_sa_sys+2.*p_sa_dias)./3;

p_pa_sys  = 0.99*p_rv_sys;
p_pa_dias = 0.27*p_rv_sys;%0.32*p_rv_sys;
p_pa_mean = (p_pa_sys+2.*p_pa_dias)./3;

p_sys_cap = 0.15*p_sa_mean; %0.26
p_pul_cap = 0.2*p_pa_mean; %0.66

p_sv_mean =  0.26*p_sa_mean;%0.16*p_sa_mean;
p_pv_mean =  0.20*p_pa_mean;%0.33*p_pa_mean;

p_la_sys  = 1.5*p_lv_dias;
p_la_dias = 0.25*p_pv_mean;

p_ra_sys  = 1.5*p_rv_dias;
p_ra_dias = 0.25*p_sv_mean;%1.*p_conv;%


%% Now, initialize the volume estimates
% These are the UNS
% Using Lumens 2009 and Boron pg. 878
Vsa_tot   = 0.14*Vtot;   % mL
Vsv_tot   = 0.70*Vtot;   % mL
Vra_tot   = 0.010*Vtot;  % mL
Vrv_tot   = 0.026*Vtot;  % mL
Vpa_tot   = 0.026*Vtot;  % mL
Vpv_tot   = 0.062*Vtot;  % mL
Vla_tot   = 0.010*Vtot;  % mL
Vlv_tot   = 0.026*Vtot;  % mL
Vsw_tot   = 0.0084*Vtot; % mL

Vm_SW = Vsw_tot;  % mL
ym    = 2.62957; % Initial guess for midwall junction (updated later)

%% Calculate stressed volumes based on Beneken and DeWitt 1966

Vsa = Vsa_tot.*0.27;
Vsv = Vsv_tot.*0.075;
Vra = Vra_tot.*1.0;
Vrv = Vrv_tot.*1.0;
Vpa = Vpa_tot.*0.58;
Vpv = Vpv_tot.*0.11;
Vla = Vla_tot.*1.0;
Vlv = Vlv_tot.*1.0;
Vsw = Vsw_tot.*1.0;

%% Next, define TriSeg parameters

Am_ref_LA = 25;%0.15;%Am_rat*Am_ref(1);
Am_ref_LV = 80;%Am_rat*Am_ref(2);
Am_ref_RA = 25;%0.15;%Am_rat*Am_ref(3);
Am_ref_RV = 100;%Am_rat*Am_ref(4);
Am_ref_SW = 35;%Am_rat*Am_ref(5);

Vwall_scale = 1.0;

V_LAwall = 1.*7.5.*Vwall_scale;%*ML2L;%7.5*ML2L;%10*ML2L;   % mL = cm^3
V_LVwall = 95.*Vwall_scale;%*ML2L;%75*ML2L;   % mL = cm^3
V_RAwall = 1.*3.5.*Vwall_scale;%*ML2L;% 3*ML2L;   % mL = cm^3
V_RVwall = 27.*Vwall_scale;%*ML2L;%30*ML2L;   % mL = cm^3
V_SWwall = 24.*Vwall_scale;%*ML2L;%31*ML2L;%40*ML2L;   % mL = cm^3

% Olsen
V0heart = (Vla+Vlv+Vra+Vrv).*1.0;


%% Parameters for the sarcomere model

% % Atria
Ls_ref  = 2.0;              % micrometer
Ls_iso  = 0.04;             % micrometer
vmax    = 12;%24;               % micrometer per second
Lsc0    = 1.51;             % micrometer
C_rest  = 0.02;             % dimensionless
tauR    = 0.03;         % seconds
tauD    = 0.03;          % seconds
tauSC   = 0.12;           % seconds
sig_act = 35; % KPa
Ls_ref_pas   = 1.8;           % micrometer
sig_pas    = 0.20;%5e-3;%0.0032; %Kpa
k_pas      = 10;%10.3207; % dimensionless

t_atr_offset = 0.15;

Spars_A = [Ls_ref; Ls_iso; vmax; Lsc0; C_rest; ...
    tauR; tauD; tauSC; sig_act; ...
    Ls_ref_pas; sig_pas; k_pas; ...
    t_atr_offset];



% Ventricles
Ls_ref  = 2.0;                % micrometer
Ls_iso  = 0.04;               % micrometer
vmax    = 12;              % micrometer per second
Lsc0    = 1.51;               % micrometer
C_rest  = 0.02;               % dimensionless
tauR    = 0.03; % seconds
tauD    = 0.08; % seconds
tauSC   = 0.40; % seconds
sig_act = 50;%75;       % KPa
Ls_ref_pas   = 1.8;           % micrometer
sig_pas    = 0.20;%5e-3;%0.0032; %Kpa
k_pas      = 10;%10.3207; % dimensionless

Spars_V = [Ls_ref; Ls_iso; vmax; Lsc0; C_rest; ...
    tauR; tauD; tauSC; sig_act; ...
    Ls_ref_pas; sig_pas; k_pas];

%% Pericardium
s = 10;%5.0; %Jezek
Peri = [V0heart; s];

%% Nominal resistance values
Ra_val = (p_lv_sys - p_sa_sys)./CO_ml;   %0.5.*p_conv/CO_ml;% KPa s / muL, Aortic Valve Resistance
Rm_val = (p_la_sys - p_lv_dias)./CO_ml;  %0.5.*p_conv/CO_ml;% KPa s / muL, Mitral Valve Resistance
Rp_val = (p_rv_sys - p_pa_sys)./CO_ml;   %0.5.*p_conv/CO_ml;% KPa s / muL, Pulmonic Valve Resistance
Rt_val = (p_ra_sys - p_rv_dias)./CO_ml;  %0.5.*p_conv/CO_ml;% KPa s / muL, Tricuspid Valve Resistance
Rvc    = (p_sv_mean - p_ra_dias)./CO_ml; % KPa s / muL, Vena Cava Resistance
Rpv    = (p_pv_mean - p_la_dias)./CO_ml; % KPa s / muL, Pulmonary venous Resistance
Rs     = (p_sa_sys-p_sys_cap)./CO_ml;   % KPa s / muL, Systemic vascular Resistance
Rp     = (p_pa_sys-p_pul_cap)./CO_ml;   % KPa s / muL, Pulmonary vascular Resistance
Csa    = Vsa./p_sa_sys;                  % muL / KPa, Systemic artery Compliance
Csv    = Vsv./p_sv_mean;                 % muL / KPa, Systemic venous Compliance
Cpa    = Vpa./p_pa_sys;                  % muL / KPa, Pulmonary artery Compliance
Cpv    = Vpv./p_pv_mean;                 % muL / KPa, Pulmonary venous Compliance

Csa = SV./(p_sa_sys-p_sys_cap);
Cpa = SV./(p_pa_sys-p_pul_cap);
% Set parameter values
Vpars = [V_LAwall; V_LVwall; V_RAwall; V_RVwall; V_SWwall; ...
    Am_ref_LA; Am_ref_LV; Am_ref_RA; Am_ref_RV; Am_ref_SW];


CVpars = [Ra_val; Rm_val; Rp_val; Rt_val; Rvc; Rpv; Rs; Rp; ...
    Csa; Csv; Cpa; Cpv];

%% Growth parameters based on Witzenburg 2018
% Eff_star = 0.17;%0.801;
% Err_star = -0.02;%-0.112;
% r_pos    = 36.4;  %r_f_pos
% r_neg    = 576;  %r_f_neg
% st50_pos = 0.097; %postst_50
% st50_neg = 0.034; %negst_50
% ff_max   = 0.1; %f_ff_max
% fr_max   = 0.1; %f_cc_max
% f_f      = 31.0; %f_f
% stim_l50 = 0.215;%sl_50
% GR = [Eff_star, Err_star, r_pos, r_neg,st50_pos, st50_neg,...
%     ff_max,fr_max,f_f,stim_l50];

%% Growth parameters from Oomen 2022 for the Hill models
Eff_star = 0.17;%0.801;
Err_star = -0.02;%-0.112;
n_pos     = 3;
n_neg     = 9;
s50_pos  = 0.075;
s50_neg  = 0.110;
f_max_pos = 0.1;
f_max_neg = 0.03;
GR = [Eff_star, Err_star, n_pos, n_neg,...
    s50_pos, s50_neg, f_max_pos,f_max_neg];
%%

if any(CVpars<0)
    error('Negative Parameters');
end

pars.V     = Vpars;
pars.SarcA = Spars_A;
pars.SarcV = Spars_V;
pars.CV    = CVpars;
pars.Peri  = Peri;
pars.GR    = GR;
pars.MI = 1; %Assume no MI at start
pars.T     = T;

%%
% To get initial values for the septal location states and the sarcomere
% states, we need to solve for zero tension to begin with.
Vm_LV = - Vlv - 0.5.*(V_LVwall+V_SWwall)+Vsw;
Vm_RV =   Vrv + 0.5.*(V_RVwall+V_SWwall)+Vsw;
Vm_LA = 0.5*V_LAwall+Vla;
Vm_RA = 0.5*V_RAwall+Vra;
% Vm = [Vm_LA Vm_LV Vm_RA Vm_RV Vm_SW];
Vm = [Vm_LV Vm_RV Vm_SW];
[xm,Am,Cm] = get_wallsegment(Vm,ym);

cla0 = C_rest;
clv0 = C_rest;
cra0 = C_rest;
crv0 = C_rest;
csw0 = C_rest;

Lla0 = Ls_ref;
Llv0 = Ls_ref;
Lra0 = Ls_ref;
Lrv0 = Ls_ref;
Lsw0 = Ls_ref;

ym0 = ym;

% Initial guess of initial conditions
IC = [Vsa; Vsv; Vra; Vrv; Vpa; Vpv; Vla; Vlv; Vsw;...
    cla0; clv0; cra0; crv0; csw0; ...
    Lla0; Llv0; Lra0; Lrv0; Lsw0; ...
    ym0];

% Update the septal midwall volume and junction point by solving a root
% finding problem
% Vm = Vm([2 4 5]);
[ym0,Vm,~,~] = get_initial_conditions(0,xm,Vm,ym,IC,pars,[Vlv;Vrv]);
IC(9)  = Vm(3);
IC(20) = ym0;

end
