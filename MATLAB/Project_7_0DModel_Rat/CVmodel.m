%%%
% this function solves the model and returns various error metrics:
% rout - residuals
% J - sum of square cost (one number)
% SS - sum of squares cost (2 diminsional vector, one entry for pressure, one for volume)
% S2 - variance estimate
% Pout - pressure model solution
% Vout - volume model solution
%%%

function [rout, J, SS, S2, N, Pout, Vout] = CVmodel(pars,data)

global ODE_TOL 

pars = exp(pars);
Init = data.Init;

%Heart parameters
Tsf = pars(8);
Trf = pars(9);
Em  = pars(10); 
EM  = pars(11);

%Initialize empty vectors for solutions 
plvS = [];  %left ventricle
VlvS = []; %left ventricle
    
T = diff(data.t_per);
NC = length(T);
tstart = 0;
tend = T(1);
i = 1;
EPS = 1e-6;

while tstart < tend
    clear plv
    clear Vlv Vao Vsa Vsv Vvc 
    clear Elv 
    
    I1 = find(data.t<tstart+EPS, 1, 'last' );
    I2 = find(data.t<tend+EPS, 1, 'last' );
    tdc = data.t(I1:I2); %current period
    
    options=odeset('RelTol',ODE_TOL, 'AbsTol',ODE_TOL);  %sets how accurate the ODE solver is
    sol = ode15s(@modelBasic,tdc,Init,options,pars,tdc(1),T(i)); %the ODE solver calls modelBasic and enters in all the following values
    sols= deval(sol,tdc); %takes the solutions at each time step in the period
    
    %assigns each row as a temporary vector to store the solutions
    %for current time period
    Vao  = sols(1,:)'; %aortic arch
    Vsa  = sols(2,:)'; %arteries
    Vsv  = sols(3,:)'; %veins
    Vvc  = sols(4,:)'; %vena cava    
    Vlv  = sols(5,:)'; %left ventricle
    
    %Determines the elasticity of the left ventricle at each period timestep
    Elv = zeros(1,length(tdc));
    for j = 1:length(tdc)
        Elv(j) = ElastanceBasic(tdc(j)-tdc(1),T(i),Tsf,Trf,Em,EM);
    end
    %Pressure of the left ventricle
    plv  = Elv'.*Vlv; 
    
    plvS  = [plvS  plv(1:end-1)']; 
    VlvS  = [VlvS   Vlv(1:end-1)'];
    
    Init = [Vao(end) Vsa(end) Vsv(end) Vvc(end) Vlv(end)];
        
    if i< NC
        tstart = tend;
        tend = tend + T(i+1);
        i = i+1;
    else
        tstart = tend;
    end
end

plvS  = [plvS   plv(end)]; 
VlvS  = [VlvS   Vlv(end)'];

%%% Error metrics
tcost = data.t(end)-.5; %time point to begin evaluating cost
ID = find(data.t_per<tcost+EPS, 1, 'last' );
tp = data.t_per(ID);
ID = find(data.t<tp+EPS, 1, 'last' );
N = length(plvS(ID:end));

VMax = mean(data.VMax);
Vmin = mean(data.Vmin);
PMax = mean(data.PMax);

% residuals normalized to magnitude of data
rp    = (plvS(ID:end) -  data.P(ID:end)')/mean(PMax);
rv    = (VlvS(ID:end) -  data.V(ID:end)')/(mean(VMax)-mean(Vmin));

% model solutions used for UQ
Pout = plvS(ID:end);
Vout = VlvS(ID:end);

SS = [rp*rp' rv*rv']; % sum of squares
J = sum(SS); %global model cost (add pressure and volume SS together)

S2 = SS/N; %variance estimate - Note, for uncertainty quanitification will need
% to rescale - multiply by magnitude of data squared

rout = [rp rv]'; %column of residuals
