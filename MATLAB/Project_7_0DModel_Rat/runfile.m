clear; clc;
%%%
% This is the primary runfile of the model - it's primary purpose is to
% solve the model and keep track of all relevant variables - including
% volume, flow, and pressure.

% Define global variables
global ODE_TOL  % Parameter describing how accurate the ODEs are solved
global Rmv Rav  % Parameters describing the values of the open and closed valuves

%load data (choose a dataset)
tshift = -0.04;  % data 1
%tshift = -0.02;  % data 2
%tshift = -0.01;  % data 3

D0 = data_process(0,1);  %data stuct w/o shifting (dataset 1, 2, 3)

% Run with nominal parameters
data = data_process(tshift,1); %make sure dataprocess#, tshift, and the .mat file loaded in line 19 all correspond to the same data set (1,2, or 3)
[x0,Init,low,hi] = load_global(data); % returns data set specifc log scaled parameters, 
pars = exp(x0); %x0 are the nominal parameters

% Run with optimized parameters
%data = data_process(tshift,3); %make sure dataprocess#, tshift, and the .mat file loaded in line 19 all correspond to the same data set (1,2, or 3)
%load('shift1/shift_-0.04_1.mat', 'x') %optimzatized parameters for data 1
%load('shift2/shift_-0.02_2.mat', 'x') %optimzatized parameters for data 1
%load('shift3/shift_-0.01_3.mat', 'x') %optimzatized parameters for data 1
%pars = exp(x);


%Resistances
Ra  = pars(1);
Rs  = pars(2); 
Rv	= pars(3);

%Elastances
Eao = pars(4);
Esa = pars(5);
Esv = pars(6);
Evc = pars(7);

%Heart parameters
Tsf  = pars(8);
Trf  = pars(9);
Emin = pars(10); 
EMax = pars(11);

%Initialize empty vectors for solutions for pressure
paoS = []; %aortic arch
psaS = []; %systemic arteries
psvS = []; %systemic Veins
pvcS = []; %Vena cava
plvS = []; %left ventricle

plvMaS = []; %max pressure per pulse 
plvmiS = []; %min pressure per pulse
plveS  = []; %end diastolic pressure per pulse

%Initialize empty vectors for solutions for volume
VaoS  = []; %aortic arch
VsaS  = []; %systemic arteries
VsvS  = []; %systemic Veins
VvcS  = []; %Vena Cava
VlvS  = []; %left ventricle
VtotS = []; %total volume

VlvMaS = []; %max volume per pulse
VlvmiS = []; %min volume per pulse
    
%Initialize empty vectors for solutions for flow
qaoS = []; %ventricle through artic arch
qvcS = []; %vena cava back to left ventricle
qaS  = []; %aortic arch to arteries
qsS  = []; %systemic flow
qvS  = []; %veins to vena cava

T = diff(data.t_per); %length of cardiac cycles
NC = length(T); %number of caridac cycles
tstart = 0;
tend = T(1);
i = 1;
EPS = 1e-6;

% this loop solves the system of ODE's iterated by each cardiac cycle -
% the system is numerically unstable if one attempts to solve the ODEs
% over multiple cycles. Within the loop volume, flow, and pressure are
% calculated and stored in appropriate variables.
while tstart < tend
    clear pao psa psv pvc plv
    clear Vlv Vao Vsa Vsv Vvc Vtot 
    clear qao qa qs qv qvc
    clear Elv 
    
    I1 = find(data.t<tstart+EPS, 1, 'last' );
    I2 = find(data.t<tend+EPS, 1, 'last' );
    tdc = data.t(I1:I2); %time vector for the current period
    %disp([tdc(1) tdc(end) T(i)]);

    options=odeset('RelTol',ODE_TOL, 'AbsTol',ODE_TOL);  %sets how accurate the ODE solver is
    sol = ode15s(@modelBasic,tdc,Init,options,pars,tdc(1),T(i)); %the ODE solver calls modelBasic and enters in all the following values
    sols= deval(sol,tdc); %deval interpolates results from ode15s back to the time vector of the data
    
    %assigns each row as a temporary vector to store the solutions
    %for current time period
    Vao  = sols(1,:)'; %aortic arch
    Vsa  = sols(2,:)'; %arteries
    Vsv  = sols(3,:)'; %veins
    Vvc  = sols(4,:)'; %vena cava    
    Vlv  = sols(5,:)'; %left ventricle
    Vtot = Vao + Vsa + Vsv + Vvc + Vlv;
    
    %Calculates volume of blood through the Pressure-Volume relationship
    %for current period
    
    pao  = Vao*Eao; %aortic arch
    psa  = Vsa*Esa; %systemic arteries
    psv  = Vsv*Esv; %systemic Veins
    pvc  = Vvc*Evc; %Vena Cava
    
    %Determines the elasticity of the left ventricle at each period timestep
    Elv = zeros(1,length(tdc));
    for j = 1:length(tdc)
        Elv(j) = ElastanceBasic(tdc(j)-tdc(1),T(i),Tsf,Trf,Emin,EMax);
    end
    %Pressure of the left ventricle
    plv  = Elv'.*Vlv; 
    
    %Flows defined by Ohm's Law
    qao  = (plv - pao)./Rav;
    qvc  = (pvc - plv)./Rmv;
    qs   = (psa  - psv)./Rs;
    qv   = (psv - pvc)./Rv;
    qa   = (pao - psa)./Ra;
    
    % Adds to pressure solution vectors every iteration of the loop
    paoS  = [paoS  pao(1:end-1)'];
    psaS  = [psaS  psa(1:end-1)'];
    psvS  = [psvS  psv(1:end-1)'];
    pvcS  = [pvcS  pvc(1:end-1)'];
    plvS  = [plvS  plv(1:end-1)']; 

    % Adds to volume solution vectors every iteration of the loop
    VaoS  = [VaoS   Vao(1:end-1)'];
    VvcS  = [VvcS   Vvc(1:end-1)'];
    VsaS  = [VsaS   Vsa(1:end-1)'];
    VsaS  = [VsvS   Vsv(1:end-1)'];
    VlvS  = [VlvS   Vlv(1:end-1)'];
    VtotS = [VtotS  Vtot(1:end-1)'];
    
    % Adds to flow solution vectors every iteration of the loop
    qaoS = [qaoS  qao(1:end-1)'];
    qvcS = [qvcS  qvc(1:end-1)'];
    qsS  = [qsS   qs(1:end-1)'];
    qvS  = [qvS   qv(1:end-1)'];
    qaS  = [qaS   qa(1:end-1)'];
    
    Init = [Vao(end) Vsa(end) Vsv(end) Vvc(end) Vlv(end)];
     
    VlvmiS(i)  = min(Vlv);
    VlvMaS(i)  = max(Vlv);
    
    plveS(i)   = plv(end);
    plvmiS(i)  = min(plv); 
    plvMaS(i)  = max(plv);
    
    if i< NC
        tstart = tend;
        tend = tend + T(i+1);
        i = i+1;
    else
        tstart = tend;
    end
end

%After for loop is done, saves last value for each of these vectors
paoS  = [paoS   pao(end)];
psaS  = [psaS   psa(end)];
psvS  = [psvS   psv(end)];
pvcS  = [pvcS   pvc(end)];
plvS  = [plvS   plv(end)]; 

VaoS  = [VaoS   Vao(end)'];
VvcS  = [VvcS   Vvc(end)'];
VsaS  = [VsaS   Vsa(end)'];
VsaS  = [VsvS   Vsv(end)'];
VlvS  = [VlvS   Vlv(end)'];
VtotS = [VtotS  Vtot(end)'];

% Adds to flow solution vectors every iteration of the loop
qaoS = [qaoS  qao(end)'];
qvcS = [qvcS  qvc(end)'];
qvS  = [qvS   qv(end)'];
qaS  = [qaS   qs(end)'];
qsS  = [qsS   qs(end)'];

figure(1);
subplot(2,1,1)
h=plot(data.t(1:10:end),plvS(1:10:end),'b',data.t(1:10:end),data.P(1:10:end),'r',D0.t(1:10:end),D0.P(1:10:end),'k:');
set(h,'Linewidth',2);
set(gca,'Fontsize',18);
ylabel('Pressure (mmHg)');
legend('model', [num2str(tshift),' shift'],'original data')
xlim([0 1])
grid on;
subplot(2,1,2)
h=plot(data.t(1:10:end),VlvS(1:10:end),'b',data.t(1:10:end),data.V(1:10:end),'r',D0.t(1:10:end),D0.V(1:10:end),'k:');
set(h,'Linewidth',2);
set(gca,'Fontsize',18);
ylabel('Volume (\muL)');
xlabel('Time (sec)')
%legend('model', [num2str(tshift),' shift'],' original data')
xlim([0 1])
grid on;

