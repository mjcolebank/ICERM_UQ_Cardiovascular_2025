% This function initializes the parameters for the model and sets initial
% values for the variables. Parmeters and initial conditions are data set
% specific.
function [x0, Init, low, hi] = load_global(data)

global DIFF_INC ODE_TOL
global Rav Rmv

ODE_TOL  = 1e-8; %ODE solving tolerance
DIFF_INC = 1e-4; %square root of the ODE tolerance - used for sensitivity analysis and optimiation

tdata_per = data.t_per;
Vdata = data.V;
Pdata = data.P;

VlvM   = mean(data.VMax);            % Left ventricular max volume FROM DATA
Vlvm   = mean(data.Vmin);            % Left ventricular min volume FROM DATA
plvM   = mean(data.PMax);            % Left ventricular max pressure FROM DATA
plvm   = mean(data.Pmin);            % Left ventriular min pressure FROM DATA

Weight   = 339;       % Weight in grams - Rat 12 - THIS IS NOT CONTAINED IN THE data STRUCT - NEEDS TO MANUALLY ENTERED HERE
HR       = round(1./mean(diff(tdata_per))) ;% beats/sec form data 
Vstr     = VlvM-Vlvm;  % ul/beat VlvM - Vlvm from data
TotalVol = 57*Weight;  % microliter [57 ml/kg = microliter/g] (weight e.g. 285g FROM DATA)
TotFlow  = HR*Vstr;    % HR * Stroke Volume 


% Flows (related to subject)
qao  = TotFlow; % flow through aortic arch to systemic arteries
qvc  = TotFlow; % flow through vena cava to heart
qa   = TotFlow; % flow through systemic arteris
qs   = TotFlow; % flow from systemic arteris to systemic veins
qv   = TotFlow; % flow from veins to vena cava


% Pressures (related to subject)     
plv    = plvM;             % Left ventricle pressure 
pao    = plv*.99;          % Aortic arch pressure
psa    = pao*.99;          % systemic arterial pressure
pvc    = plvm*1.1;         % Pressure in Vena Cave
psv    = pvc *1.1;         % Systemic Vein Pressure
paod   = pao-30;           % diasoltic aortic pressure
psad   = psa-30;           % diasoltic systemic artery pressure
paom   = paod + 1/3*(pao-paod); % mean aortic pressure
psam   = psad + 1/3*(psa-psad); % mean systemic artery pressure

% Resistances (Ohm's law)
Rs  = (psam - psv)/qs;     
Ra  = (paom - psam)/qa;
Rv  = (psv - pvc)/qv;
Rav = (plvM - pao)/qao; 
Rmv = (pvc - plvm)/qvc;

% Volumes (Beneken and deWit)
Vsa  = TotalVol*0.20;   
Vsv  = TotalVol*0.70;   
Vvc  = TotalVol*0.075;
Vao  = TotalVol*0.025;

%stressed volumes (systole)
VsaS = Vsa*0.3;
VsvS = Vsv*0.3;   
VvcS = Vvc*0.3;
VaoS = Vao*0.3;

% Elastances (Beneken)
Eao = pao/VaoS;      % 
Esa = psa/VsaS;      % 
Esv = psv/VsvS;      % 
Evc = pvc/VvcS;      % 

% Ventricle parameters 
% Average of all segments - these values were determined by inspection
Tsf = 0.25;     % Time fraction for systolic phase, increasing elasticity 
Trf = 0.60;     % Time fraction for Systolic phase, decreasing elasticity 

Emin= pvc/VlvM;  % Minimum elasticity of the heart [ plvm/(VlvM-Vd) ]
Emax= paod/Vlvm; % Maximum elasticity of the heart [ pao_diastole/(Vlvm-Vd) ] Pressure diastole estimated

% Initial conditions
Init = [VaoS VsaS VsvS VvcS VlvM];

% Parameter vector
x0 = [Ra Rs Rv...
      Eao Esa Esv Evc...
      Tsf Trf Emin Emax]';   %8-11
x0 = log(x0);

% Upper bound and lower bounds for optimization
hi  = x0+log(20);   % x0*20 (when not log)
low = x0-log(20);   % x0/20 (when not log)


