function xdot = modelBasic(t,y,pars,Vunstr,tst,T)

Vsa    = y(1);  % systemic artery volume
Vsv    = y(2);  % systemic venous volume
Vpa    = y(3);  % pulmonary artery volume
Vpv    = y(4);  % pulmonary venous volume
Vra    = y(5);  % right atria volume
Vla    = y(6);  % left atria volume
Vrv    = y(7);  % right ventricular volume
Vlv    = y(8);  % left ventricular volume

% Unstressed Volumes
VsaU = Vunstr(1);
VsvU = Vunstr(2);
VpaU = Vunstr(3);
VpvU = Vunstr(4);
VaU  = Vunstr(5);
VvU  = Vunstr(6);

%% PARAMETERS
% Resistances
Rs    = pars(1); % systemic periphery
Rp    = pars(2); % pulmonary periphery
Rava  = pars(3); % aortic valve
Rmva  = pars(4); % mitral valve
Rpva  = pars(5); % pulmonary valve
Rtva  = pars(6); % tricuspid valve
Rpv   = pars(7); % pulmonary veins
Rsv   = pars(8); % systemic veins

% Compliances
Csa   = pars(9); % systemic arteries
Csv   = pars(10); % systemic veins
Cpa   = pars(11); % pulmonary arteries
Cpv   = pars(12); % pulmonary veins

% Heart
% Elastance
EMra  = pars(13); % right atrium max 
Emra  = pars(14); % right atrium min 
EMla  = pars(15); % left atrium max 
Emla  = pars(16); % left atrium min 
EMrv  = pars(17); % right ventricle max 
Emrv  = pars(18); % right ventricle min 
EMlv  = pars(19); % left ventricle max
Emlv  = pars(20); % left ventricle min 

% Timing
Trra  = pars(21); % end right atrium relaxation
tcra  = pars(22); % right atrium begin contraction
Tcra  = pars(23); % right atrium end contraction

Trla  = Trra;%*1.01;  % end left atrium relaxation
tcla  = tcra;%*1.05;  % left atrium begin contraction
Tcla  = Tcra;       % left atrium end contraction

Tcrv  = pars(24);    % right ventricle contraction
Trrv  = pars(25);    % right ventricle relaxation
Tclv  = Tcrv;%*0.95;   % left ventricle contraction
Trlv  = Trrv;        % left ventricle relaxation

% Timevarying elastance
Era  = ElastanceAtrium(t-tst,EMra,Emra,Trra,tcra,Tcra,T); % Elatance right atrium
Ela  = ElastanceAtrium(t-tst,EMla,Emla,Trla,tcla,Tcla,T); % Elastance left atrium
Erv  = ElastanceVentricle(t-tst,EMrv,Emrv,Tcrv,Trrv);     % Elastance right ventricle
Elv  = ElastanceVentricle(t-tst,EMlv,Emlv,Tclv,Trlv);     % Elastance left ventricle

% Pressuref
psa  = (Vsa-VsaU)/Csa; % systemic arteries 
psv  = (Vsv-VsvU)/Csv; % systemic veins 
ppa  = (Vpa-VpaU)/Cpa; % pressure pulmonary artery
ppv  = (Vpv-VpvU)/Cpv; % pressure pulmonary vein

% % Using new activation function
plv  = Elv*(Vlv-VvU); % pressure left ventricle
prv  = Erv*(Vrv-VvU); % pressure right ventricle
pla  = Ela*(Vla-VaU); % pressure left ventricle
pra  = Era*(Vra-VaU); % pressure right ventricle
         
if plv > psa % flow through aortic valve 
    qava = (plv-psa)/Rava;
else
    qava  = 0;
end

if pla > plv % flow through mitral valve
    qmva = (pla-plv)/Rmva; % valve open
else
    qmva = 0;             % valve closed
end

if prv > ppa % flow through pulmonary valve
    qpva = (prv-ppa)/Rpva;
else
    qpva  = 0;
end

if pra > prv  % flow through tricuspid valve
    qtva = (pra-prv)/Rtva; % valve open
else
    qtva = 0;             % valve closed
end

if psv > pra  % flow through tricuspid valve
    qsv = (psv-pra)/Rsv;
else
    qsv = 0;
end

% Remaining flows
qpv  = (ppv - pla) / Rpv;
qs   = (psa - psv) / Rs; % Systemic periphery
qp   = (ppa - ppv) / Rp; % Pulmonary periphery


% Differential Equations
dVsa = qava - qs;
dVsv = qs   - qsv;
dVpa = qpva - qp;
dVpv = qp   - qpv;
dVra = qsv  - qtva; 
dVla = qpv  - qmva;
dVrv = qtva - qpva;
dVlv = qmva - qava;

xdot = [dVsa;...
        dVsv;...
        dVpa;...
        dVpv;...
        dVra;...
        dVla;...
        dVrv;...
        dVlv;];
end
