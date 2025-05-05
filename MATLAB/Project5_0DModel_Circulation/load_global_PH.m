% This function initializes the parameters for the model and sets initial
% values for the variables.
function [pars,Vunstr,Init,low,hi,data] = load_global_PH(data,ptnb)

% Total blood volume (mL)
BSA = sqrt(data.Height*data.Weight/3600);
try 
  BSA = data.BSA;
catch
end

if data.Gender == 2
  TotalVol = (3.47*BSA - 1.954)*1000; % female
elseif data.Gender == 1
  TotalVol = (3.29*BSA - 1.229)*1000; % male
end

% Cardiac Output (mL/s)
qtot = data.CO*1000/60;

% Stroke Volume
Vstroke = qtot*data.T;   % (mL/s * s/beat) = mL/beat)
data.Vstr = Vstroke;

% Blood Pressure
% Systemic arteries
psa    = data.pSAa;      % mean systemic artery pressure
pSsa   = data.pSAM;      % systolic systemic artery pressure
pDsa   = data.pSAm;      % systolic systemic artery pressure

% Pulmonary arteries
ppa    = mean(data.pPA); % mean pulmonary artery pressure (area under curve) 
pSpa   = max(data.pPA);  % systolic pulmonary artery pressure
pDpa   = min(data.pPA);  % systolic pulmonary artery pressure
pmra   = min(data.pRA(1:round(length(data.pRA)/3)));

% Veins
ppv    = min(data.pPW);  % pulmonary veins
psv    = 1.10*pmra;      % systemic veins

% The heart
% Atria
Nra = round(2*length(data.pRA)/3);
pMra   = max(data.pRA(Nra:end));  % max right atrial pressure
pmla   = 0.9*ppv;                 % min left atrial pressure
pMla   = pmla + 5;                % max left atrial pressure   

% Ventricles
pMrv   = max(data.pRV); % max right ventricular pressure
pmrv   = min(data.pRV); % min right ventricular pressure
pMlv   = 1.02*pSsa;     % max left ventricular pressure
pmlv   = 0.98*pmla;     % min left ventricular pressure

% Blood Volumes (From Tello et al. 2019,doi:10.1152/ajpheart.00485.2019)
VraM  = BSA*(58.9);     % max right atrial volume
Vram  = BSA*(36.2);     % min right atrial volume
VlaM  = BSA*30;         % max left atrial volume
Vlam  = BSA*13;         % min left atrial volume

VrvM = BSA*(116.9);     % max right ventricular volume
VlvM = BSA*(80);        % max left ventricular volume
Vrvm = VrvM-Vstroke;    % min right ventricular volume
Vlvm = VlvM-Vstroke;    % min left ventricular volume

VraMF = VraM/TotalVol;  % max right atrial volume fraction
VlaMF = VlaM/TotalVol;  % max left atrial volume fraction
VrvMF = VrvM/TotalVol;  % max right ventricle volume fraction
VlvMF = VlvM/TotalVol;  % max left ventricle volume fraction
HVF   = VraMF+VlaMF+VrvMF+VlvMF; % heart volume fraction

% Volume distribution
Vsa = 0.13*TotalVol;       % Arterial volume (13%)
Vsv = (0.72-HVF)*TotalVol; % Venous volume (~65%)  
Vpa = 0.03*TotalVol;       % Pulmonary artery (3%)
Vpv = 0.11*TotalVol;       % Pulmonary vein volume (11%)
 
% Stressed Volumes (from Beneken)
VsaU = Vsa*(1-0.27);   % Unstressed systemic arteries volume 73% 
VsvU = Vsv*(1-0.08);   % Unstressed systemic venous volume 92% 
VpaU = Vpa*(1-0.58);   % Unstressed pulmonary arterial volume 62%
VpvU = Vpv*(1-0.11);   % Unstreseed pulmonary venous volume 89%
VaU  = 10;             % Unstressed atrial volume 10mL
VvU  = 5;              % Unstressed ventricular volume 5 mL
Vunstr = [VsaU VsvU VpaU VpvU VaU VvU]; % Vector of unstressed volumes

% Volume vector
data.VlvM = VlvM; % max left ventricular volume
data.VrvM = VrvM; % max right ventricular volume
data.VlaM = VlaM; % max left atrial volume
data.VraM = VraM; % max right atrial volume
data.Vlam = Vlam; % min left atrial volume
data.Vram = Vram; % min right atrial volume

data.RaVstr = VraM-Vram; % right atrial stroke volume
data.LaVstr = VlaM-Vlam; % left atrial stroke volume

% Resistence
% Peripheral resistances (Ohm's law)
Rs = (pDsa-psv)/qtot;    % systemic resistance 
Rp = (pDpa-ppv)/qtot;    % pulmonary resistance

% Valve resistances
Rava = (pMlv-pSsa)/qtot; % resistance aortic valve
Rpva = (pMrv-pSpa)/qtot; % resistance pulmonary valve 
Rmva = 0.01;             % resistance mitral valve
Rtva = 0.03;             % resistance tricuspid valve

% Venous resistances
Rpv = (ppv-pmla)/qtot;   % resistance of pulmonary vein
Rsv = (psv-pmra)/qtot;   % resistance of vena cava

% Compliance
Csa = (Vsa-VsaU)/pDsa;   % systemic artery compliance
Csv = (Vsv-VsvU)/psv;    % systemic venous compliance  
Cpa = (Vpa-VpaU)/pDpa;   % pulmonary artery compliance
Cpv = (Vpv-VpvU)/ppv;    % pulmonary vein compliance

% Heart parameters
% Elastance parameters atria 
Emra = pmra/(VraM-VaU);  % min right atrium 
EMra = pMra/(Vram-VaU);  % max right atrium 
Emla = pmla/(VlaM-VaU);  % min left atrium
EMla = pMla/(Vlam-VaU);  % max left atrium

% Elastance parameters ventricles
Emrv = pmrv/(VrvM-VvU);  % min right ventricle
EMrv = pMrv/(Vrvm-VvU);  % max right ventricle    
Emlv = pmlv/(VlvM-VvU);  % min left ventricle
EMlv = pMlv/(Vlvm-VvU);  % max left ventricle

% Timing heart contraction 
% Aria
Trra = data.Trra;
tcra = data.tcra;
Tcra = data.Tcra;

% Ventricle
Tcrv = data.Tcrv;
Trrv = data.Trrv;

% Initial conditions
Init = [Vsa Vsv Vpa Vpv VraM VlaM VrvM VlvM];

% Parameter vector
x0 = [Rs Rp Rava Rmva Rpva Rtva Rpv Rsv ...  % 1-8
      Csa Csv Cpa Cpv ...                    % 9-12
      EMra Emra EMla Emla ...                % 13-16
      EMrv Emrv EMlv Emlv ...                % 17-20
      Trra tcra Tcra ...                     % 21-23
      Tcrv Trrv]';                           % 24-25
pars = log(x0);

% Upper and lower bounds for optimization
hi  = pars + log(4);
low = pars - log(4);

% Prevent overlap of bounds for timing parameters
low(21) = log(x0(21)-0.1);
hi (21) = log(x0(21)+0.15);
low(22) = log(x0(22)-0.15);
hi (22) = log(x0(22)+0.15);
low(23) = log(x0(23)-0.15);
hi (23) = log(data.T);
low(24) = log(x0(24)-0.15);
hi (24) = log(x0(24)+0.15);
low(25) = log(x0(25)-0.15);
hi(25)  = log(x0(25)+0.15);

if ptnb ==1
    hi (22) = log(x0(22)+0.065);
    low(23) = log(x0(23)-0.065);
elseif ptnb==2
    hi (22) = log(x0(22)+0.054);
    low(23) = log(x0(23)-0.054);
    hi (24) = log(x0(24)+0.095);
    low(25) = log(x0(25)-0.095);
elseif ptnb==3
    hi (22) = log(x0(22)+0.09);
    low(23) = log(x0(23)-0.09);
elseif ptnb==4
    hi (22) = log(x0(22)+0.1);
    low(23) = log(x0(23)-0.1);
elseif ptnb==5
    hi (22) = log(x0(22)+0.14);
    low(23) = log(x0(23)-0.14);
    hi (24) = log(x0(24)+0.1);
    low(25) = log(x0(25)-0.1);
elseif ptnb==6
    hi (22) = log(x0(22)+0.09);
    low(23) = log(x0(23)-0.09);
    hi (24) = log(x0(24)+0.12);
    low(25) = log(x0(25)-0.12);
elseif ptnb==7
    hi (22) = log(x0(22)+0.13);
    low(23) = log(x0(23)-0.13);
elseif ptnb==8
    hi (22) = log(x0(22)+0.06);
    low(23) = log(x0(23)-0.06);
    hi (24) = log(x0(24)+0.12);
    low(25) = log(x0(25)-0.12);
elseif ptnb==9
    hi (22) = log(x0(22)+0.09); 
    low(23) = log(x0(23)-0.09); 
    hi (24) = log(x0(24)+0.12);
    low(25) = log(x0(25)-0.12);
end;

