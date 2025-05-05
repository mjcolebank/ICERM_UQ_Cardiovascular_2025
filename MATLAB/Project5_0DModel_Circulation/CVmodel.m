% solves the model and plots solutions against data
%
% This function takes in the data structure, and then runs the model until
% convergence has been met for all the model states.
function [rout,J,CV] = CVmodel(pars,data)

ODE_TOL = data.ODE_TOL;
REL_TOL = data.REL_TOL;

VUsa = data.V(1);
VUsv = data.V(2);
VUpa = data.V(3);
VUpv = data.V(4);
VaU  = data.V(5);
VvU  = data.V(6);
Vunstr = data.V;

pars = exp(pars);

Init = data.Init;

res_flag = data.res_flag;

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
Csa = pars(9); % systemic arteries
Csv = pars(10); % systemic veins
Cpa = pars(11); % pulmonary arteries
Cpv = pars(12); % pulmonary veins

% Elatance Heart parameters
EMra   = pars(13); % max right atrial elastance
Emra   = pars(14); % min right atrial elastance
EMla   = pars(15); % max left ventricle elastance
Emla   = pars(16); % min right ventricle elastance
EMrv   = pars(17); % max left ventricle elatance
Emrv   = pars(18); % min left ventricle elatance
EMlv   = pars(19); % max left ventricle elatance
Emlv   = pars(20); % min left ventricle elatance

% Heart timing parameters
% Atrium 
Trra  = pars(21); % end right atrium relaxation
tcra  = pars(22); % right atrium begin contraction
Tcra  = pars(23); % right atrium end contraction

Trla  = Trra;     % end left atrium relaxation
tcla  = tcra;     % left atrium begin contraction
Tcla  = Tcra;     % left atrium end contraction

% Ventricle
Tcrv  = pars(24);   % right ventricle contraction
Trrv  = pars(25);   % right ventricle relaxation
Tclv  = Tcrv;       % left ventricle contraction
Trlv  = Trrv;       % left ventricle relaxation

% Initialize vectors for pressure solutions
CV.plaS   = [];
CV.plvS   = [];
CV.psvS   = [];
CV.ppaS   = [];  
CV.pPAdS  = [];
CV.ppaMS  = [];
CV.pPAMdS = [];    
CV.ppamS  = [];
CV.pPAmdS = [];
CV.praS   = [];
CV.pRAdS  = [];
CV.praaMS = [];
CV.pRAaMdS= [];
CV.praamS = [];
CV.pRAamdS= [];
CV.praMS  = [];
CV.pRAMdS = [];
CV.pramS  = [];
CV.pRAmdS = [];
CV.pravMS = [];
CV.pRAvMdS= [];
CV.pravmS = [];
CV.pRAvmdS= [];
CV.prvS   = [];
CV.pRVdS  = [];
CV.prvMS  = [];
CV.pRVMdS = [];
CV.prvmS  = [];
CV.pRVmdS = [];
CV.psaS   = [];
CV.psaMS  = [];
CV.pSAMdS = [];
CV.psamS  = [];
CV.pSAmdS = [];
CV.ppvS   = [];
CV.ppvmS  = [];
CV.pPWdS  = [];
CV.COS    = [];
CV.COdS   = [];
CV.tdS    = [];

CV.VsaS   = [];
CV.VsvS   = [];
CV.VpaS   = [];
CV.VpvS   = [];
CV.VraS   = [];
CV.VrvS   = [];
CV.VlaS   = [];
CV.VlvS   = [];

CV.VlvMS  = [];
CV.VrvMS  = [];
CV.VlaMS  = [];
CV.VraMS  = [];

CV.VlamS  = [];
CV.VramS  = [];

CV.VladMS  = [];
CV.VradMS  = [];
CV.VladmS  = [];
CV.VradmS  = [];
CV.VlvdMS  = [];
CV.VrvdMS  = [];

CV.RVstrS = [];
CV.LVstrS = [];
CV.VstrdS = [];

CV.RaVstrS = [];
CV.LaVstrS = [];
CV.RaVstrdS = [];    
CV.LaVstrdS = [];    
   
CV.tperS  = [];
CV.tdS    = [];

CV.qava = [];
CV.qmva = [];
CV.qtva = [];
CV.qsv  = [];
CV.qpva = [];
    
T  = data.T;
td = data.td;
dt = data.dt;
NC = data.NC; 

M2   = round(NC/2);

k1 = 1; % index of first time step in first period
k2 = round(T/dt)+k1; %index of last time step in first period
for i = 1:NC % go through loop NC times
        
    clear Vsa Vsv Vpa Vpv  Vra Vla Vrv Vlv
    clear psa psv ppa ppv  Era Era Ela pla Erv prv Elv plv
    clear qs qp  qmva qtva qsv qpv
        
    tdc = td(k1:k2);  %current time
  
    options = odeset('RelTol',REL_TOL, 'AbsTol',ODE_TOL);       %sets how accurate the ODE solver is
    sol     = ode15s(@modelBasic,tdc,Init,options,pars,Vunstr,tdc(1),T); %solves the ODE uses modelBasic
    sols    = deval(sol,tdc);                                     %interpolates solution at data times
        
    % extract solutions
    Vsa  = sols(1,:)';  % systemic arteries
    Vsv  = sols(2,:)';  % systemic veins
    Vpa  = sols(3,:)';  % pulmonary arteries
    Vpv  = sols(4,:)';  % pulmonary veins
    Vra  = sols(5,:)';  % right atrium
    Vla  = sols(6,:)';  % left atrium
    Vrv  = sols(7,:)';  % right ventricle
    Vlv  = sols(8,:)';  % left ventricle
        
    psa = (Vsa-VUsa)/Csa; % systemic arteries
    psv = (Vsv-VUsv)/Csv; % systemic veins
    ppa = (Vpa-VUpa)/Cpa; % pressure pulmonary artery
    ppv = (Vpv-VUpv)/Cpv; % pressure pulmonary vein
        
    % Heart elastance
    
    for w = 1:length(tdc)
        Era(w)  = ElastanceAtrium(tdc(w)-tdc(1),EMra,Emra,Trra,tcra,Tcra,T); % Right atrium elastance
        Ela(w)  = ElastanceAtrium(tdc(w)-tdc(1),EMla,Emla,Trla,tcla,Tcla,T); % Left atrium elastance
        
        Erv(w)  = ElastanceVentricle(tdc(w)-tdc(1),EMrv,Emrv,Tcrv,Trrv);   % Right ventricle elastance
        Elv(w)  = ElastanceVentricle(tdc(w)-tdc(1),EMlv,Emlv,Tclv,Trlv);   % Left ventricle elastance
    end

    prv = Erv'.*(Vrv-VvU); % Right ventricle pressure
    plv = Elv'.*(Vlv-VvU); % Left ventricle pressure
    pra = Era'.*(Vra-VaU); % Right ventricle pressure
    pla = Ela'.*(Vla-VaU); % Right ventricle pressure
     
    % Flows defined by Ohm's Law
    qs    = (psa - psv)/Rs; % systemic periphery
    qp    = (ppa - ppv)/Rp; % pulmonary periiphery
    
    % Cardiac Output
    CO = trapz(tdc, qs)/(tdc(end)-tdc(1))/1000*60;
    
    % Flows
    qava = max((plv-psa)/Rava,0);
    qmva = max((pla-plv)/Rmva,0);
    qtva = max((pra-prv)/Rtva,0);
    qsv  = max((psv-pra)/Rsv,0);
    qpva = max((prv-ppa)/Rpva,0);
    qpv  = (ppv-pla)/Rpv;
    
    % time
    CV.tdS    = [CV.tdS tdc(1:end-1)];
    CV.tperS  = [CV.tperS tdc(1)];
 
    % pla
    CV.plaS   = [CV.plaS pla(1:end-1)'];
    
    % plv
    CV.plvS   = [CV.plvS plv(1:end-1)'];
    
    % psv
    CV.psvS   = [CV.psvS psv(1:end-1)'];
    
    % ppa
    CV.ppaS   = [CV.ppaS ppa(1:end-1)'];
    CV.pPAdS  = [CV.pPAdS data.pPA(1:end-1)'];
    CV.ppaMS  = [CV.ppaMS  max(ppa)];  
    CV.pPAMdS = [CV.pPAMdS max(data.pPA)];    
    CV.ppamS  = [CV.ppamS  min(ppa)];
    CV.pPAmdS = [CV.pPAmdS min(data.pPA)];
    
    % pra
    CV.praS   = [CV.praS pra(1:end-1)'];
    CV.pRAdS  = [CV.pRAdS data.pRA(1:end-1)'];
    Nra = round(2*length(pra)/3);    
    CV.praMS  = [CV.praMS  max(pra(Nra:end))];
    CV.pRAMdS = [CV.pRAMdS max(data.pRA(Nra:end))];
    CV.pramS  = [CV.pramS  min(pra)];
    CV.pRAmdS = [CV.pRAmdS min(data.pRA(1:round(length(data.pRA)/3)))];
       
    % prv
    CV.prvS   = [CV.prvS prv(1:end-1)'];
    CV.pRVdS  = [CV.pRVdS data.pRV(1:end-1)'];
    CV.prvMS  = [CV.prvMS  max(prv)];
    CV.pRVMdS = [CV.pRVMdS max(data.pRV)];
    CV.prvmS  = [CV.prvmS  min(prv)];
    CV.pRVmdS = [CV.pRVmdS min(data.pRV)];
    
    % psa
    CV.psaS   = [CV.psaS psa(1:end-1)'];
    CV.psaMS  = [CV.psaMS  max(psa)];
    CV.pSAMdS = [CV.pSAMdS data.pSAM];
    CV.psamS  = [CV.psamS  min(psa)];
    CV.pSAmdS = [CV.pSAmdS data.pSAm];

    % ppv
    CV.ppvmS  = [CV.ppvmS   mean(ppv)];
    CV.ppvS   = [CV.ppvS   ppv(1:end-1)'];
    CV.pPWdS  = [CV.pPWdS  min(data.pPW)];
    
    % CO
    CV.COS    = [CV.COS CO];
    CV.COdS   = [CV.COdS data.CO]; 
    
    % Volumes
    CV.VsaS   = [CV.VsaS Vsa(1:end-1)'];
    CV.VsvS   = [CV.VsvS Vsv(1:end-1)'];
    CV.VpaS   = [CV.VpaS Vpa(1:end-1)'];
    CV.VpvS   = [CV.VpvS Vpv(1:end-1)'];
  
    CV.VraS   = [CV.VraS Vra(1:end-1)'];
    CV.VrvS   = [CV.VrvS Vrv(1:end-1)'];
    CV.VlaS   = [CV.VlaS Vla(1:end-1)'];
    CV.VlvS   = [CV.VlvS Vlv(1:end-1)'];
    
    CV.VlvMS  = [CV.VlvMS max(Vlv)];
    CV.VrvMS  = [CV.VrvMS max(Vrv)];
    CV.VlaMS  = [CV.VlaMS max(Vla)];
    CV.VraMS  = [CV.VraMS max(Vra)];
  
    CV.VlamS  = [CV.VlamS min(Vla)];
    CV.VramS  = [CV.VramS min(Vra)];
  
    CV.VlvdMS = [CV.VlvdMS data.VlvM];
    CV.VrvdMS = [CV.VrvdMS data.VrvM];
    CV.VladMS = [CV.VladMS data.VlaM];
    CV.VradMS = [CV.VradMS data.VraM];
    
    CV.RVstrS = [CV.RVstrS max(Vrv)-min(Vrv)];
    CV.LVstrS = [CV.LVstrS max(Vlv)-min(Vlv)];
    CV.VstrdS = [CV.VstrdS data.Vstr];
    
    CV.RaVstrS = [CV.RaVstrS max(Vra)-min(Vra)];
    CV.LaVstrS = [CV.LaVstrS max(Vla)-min(Vla)];
    CV.RaVstrdS = [CV.RaVstrdS data.RaVstr];    
    CV.LaVstrdS = [CV.LaVstrdS data.LaVstr];    
    
    CV.qava = [CV.qava qava(1:end-1)'];
    CV.qmva = [CV.qmva qmva(1:end-1)'];
    CV.qtva = [CV.qtva qtva(1:end-1)'];
    CV.qsv  = [CV.qsv  qsv(1:end-1)'];
    CV.qpva = [CV.qpva qpva(1:end-1)'];

    
    Init = [Vsa(end) Vsv(end) Vpa(end) Vpv(end) Vra(end) Vla(end) Vrv(end) Vlv(end)];
    if i < NC
        k1 = k2; %sets last index of this loop as first index for next loop
        k2 = round(T/dt)+k1;
    end
    
end

I1 = max(find(tdc-tdc(1)<Trra));
I2 = max(find(tdc-tdc(1)<tcra));
I3 = max(find(tdc-tdc(1)<tcla));
    
% save solution
CV.plaS   = [CV.plaS pla(end)];
CV.plvS   = [CV.plvS plv(end)];
CV.psvS   = [CV.psvS psv(end)];
CV.ppaS   = [CV.ppaS ppa(end)];
CV.pPAdS  = [CV.pPAdS data.pPA(end)];
CV.ppvS   = [CV.ppvS ppv(end)];
CV.praS   = [CV.praS pra(end)];
CV.pRAdS  = [CV.pRAdS data.pRA(end)];
CV.prvS   = [CV.prvS prv(end)];
CV.pRVdS  = [CV.pRVdS data.pRV(end)];
CV.psaS   = [CV.psaS psa(end)];

CV.VsaS   = [CV.VsaS Vsa(end)];
CV.VsvS   = [CV.VsvS Vsv(end)];
CV.VpaS   = [CV.VpaS Vpa(end)];
CV.VpvS   = [CV.VpvS Vpv(end)];
CV.VraS   = [CV.VraS Vra(end)];
CV.VrvS   = [CV.VrvS Vrv(end)];
CV.VlaS   = [CV.VlaS Vla(end)];
CV.VlvS   = [CV.VlvS Vlv(end)];

CV.qava = [CV.qava qava(end)];
CV.qmva = [CV.qmva qmva(end)];
CV.qtva = [CV.qtva qtva(end)];
CV.qsv  = [CV.qsv  qsv(end)];
CV.qpva = [CV.qpva qpva(end)];

CV.tdS    = [CV.tdS tdc(end)]; 

M2   = round(NC/2);
N2   = round(length(CV.ppaS)/2); 

w = 1;

% MJC this is the correct penalty function
penalty_VlvM = w*max(0,(CV.VlvdMS(M2:end)-CV.VlvMS(M2:end))./CV.VlvdMS(M2:end)/sqrt(M2));
penalty_VrvM = w*max(0,(CV.VrvdMS(M2:end)-CV.VrvMS(M2:end))./CV.VrvdMS(M2:end)/sqrt(M2));
penalty_VlaM = w*max(0,(CV.VladMS(M2:end)-CV.VlaMS(M2:end))./CV.VladMS(M2:end)/sqrt(M2));
penalty_VraM = w*max(0,(CV.VradMS(M2:end)-CV.VraMS(M2:end))./CV.VradMS(M2:end)/sqrt(M2));


routS = [(CV.ppaMS(M2:end)-CV.pPAMdS(M2:end))./CV.pPAMdS(M2:end)/sqrt(M2) ...
         (CV.ppamS(M2:end)-CV.pPAmdS(M2:end))./CV.pPAmdS(M2:end)/sqrt(M2) ...
         (CV.praMS(M2:end)-CV.pRAMdS(M2:end))./CV.pRAMdS(M2:end)/sqrt(M2) ... 
         (CV.pramS(M2:end)-CV.pRAmdS(M2:end))./CV.pRAmdS(M2:end)/sqrt(M2) ...
         (CV.prvMS(M2:end)-CV.pRVMdS(M2:end))./CV.pRVMdS(M2:end)/sqrt(M2) ... 
         (CV.prvmS(M2:end)-CV.pRVmdS(M2:end))./CV.pRVmdS(M2:end)/sqrt(M2) ... 
         (CV.psaMS(M2:end)-CV.pSAMdS(M2:end))./CV.pSAMdS(M2:end)/sqrt(M2) ... 
         (CV.psamS(M2:end)-CV.pSAmdS(M2:end))./CV.pSAmdS(M2:end)/sqrt(M2) ...
         (CV.ppvmS(M2:end)-CV.pPWdS(M2:end)) ./CV.pPWdS(M2:end) /sqrt(M2) ... 
         (CV.COS(M2:end)  -CV.COdS(M2:end))  ./CV.COdS(M2:end)  /sqrt(M2) ...
         (CV.RVstrS(M2:end)-CV.VstrdS(M2:end))./CV.VstrdS(M2:end)/sqrt(M2) ...
         (CV.LVstrS(M2:end)-CV.VstrdS(M2:end))./CV.VstrdS(M2:end)/sqrt(M2) ... 
         penalty_VlvM penalty_VrvM penalty_VlaM penalty_VraM];
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
LS = 1;

if res_flag == 1
    % Static only
    routT = [];
    LT = 1;
elseif res_flag == 2
    % Static and dynamic
    routT   = [(CV.praS(N2:end)-CV.pRAdS(N2:end))./max(CV.pRAdS(N2:end))/sqrt(N2) ...
               (CV.prvS(N2:end)-CV.pRVdS(N2:end))./max(CV.pRVdS(N2:end))/sqrt(N2) ...
               (CV.ppaS(N2:end)-CV.pPAdS(N2:end))./max(CV.pPAdS(N2:end))/sqrt(N2)]; 
   LT = 1;
end;

rout = [routS/LS routT/LT]';
J    = rout'*rout;

end
