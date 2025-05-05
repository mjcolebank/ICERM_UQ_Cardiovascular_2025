function RunSimulation(R,INDMAP,task,par,ptnb,iter,rnum) 
       
%cd PatientWaveforms
cd PatientPhysiol
if ptnb == 1 
    data = patient1;
elseif ptnb == 2
   data = patient2; 
elseif ptnb == 3
   data = patient3;   
elseif ptnb == 4
   data = patient4; 
elseif ptnb == 5
   data = patient5; 
elseif ptnb == 6
   data = patient6; 
elseif ptnb == 7
   data = patient8; 
elseif ptnb == 8
   data = patient11;  
else
   data = patient20;
end
cd ..

data.INDMAP = INDMAP;

data.res_flag = R; 

data.ODE_TOL  = 1e-8;
data.REL_TOL  = 1e-8;
data.DIFF_INC = 1e-4;

[pars,Vunstr,Init,low,hi,data] = load_global_PH(data,ptnb);
pars = exp(pars);
pars = pars.*(1+0.2*rnum);
pars = log(pars);

data.Init= Init;
data.V   = Vunstr;
data.hi  = hi;
data.low = low;

if task == 1 % SENSITIVTY ANALYIS
   history = DriverBasic_sens(data,pars);
   s = strcat('Sens_',num2str(ptnb),'_R',num2str(R),'.mat');
   save(s,'history');
elseif task == 2  % OPTIMIZATION Pick appriate string
   [history, optpars] = DriverBasic_opt(data,pars);

   s = strcat('Opt_',num2str(ptnb),'_R',num2str(R),'_iter',num2str(iter),'.mat');
   save(s,'history','optpars','pars','data','R');
elseif task == 3    
  if par == 2
    s  = strcat('Opt_',num2str(ptnb),'_R',num2str(R),'_iter',num2str(iter));
    sp = strcat('Opt_',num2str(ptnb),'_R',num2str(R));
    load(strcat(s,'.mat'));
    orgp = exp(pars);
    pars(history.INDMAP) = optpars(history.INDMAP); 
    optp = exp(pars);
  end
  sp = strcat('Nom_',num2str(ptnb),'_R',num2str(R));
  DriverBasic_plot(pars,data,ptnb,sp);
elseif task == 4
    s = strcat('SensVstr2_p',num2str(ptnb),'_R',num2str(R)','.mat');
    load(s);
    plotSens(history);
end
