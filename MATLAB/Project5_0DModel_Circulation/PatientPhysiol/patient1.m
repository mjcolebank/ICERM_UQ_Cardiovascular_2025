    
% Data file for patient 1:
% This file contains all the information needed to run any of the
% driverbasic functions for the given patient.

function data = patient1

% Load the original data file
load('Patient_1_waveforms.mat');
load('Patient_1.mat');
load('Patient_1_newW.mat');

HR=[69 68 68 60 68 65 65 66 64 64 65];
HRa = mean(HR);
T = 60/HRa;

% Reassign the timing parameters
dt  = 0.01;
T   = round(T/dt)*dt;
td_per = 0:dt:T;

NC  = 30; 
td  = 0:dt:T*NC;

% end and beginning pressures should be the same
pRA = [pRA pRA(1)];
pRV = [pRV pRV(1)];
pPW = [pPW pPW(1)];
pPA = [pPA pPA(1)];

tdata   = [0:4096]/4097*T;
tdataRA = [0:3563]/3564*T;

% Shift RA
pRA = circshift(pRA,-950);

% Smoothing and interpolating
pPAi = interp1(tdata,pPA,td_per);
pPWi = interp1(tdata,pPW,td_per);
pRAi = interp1(tdataRA,pRA,td_per);
pRVi = interp1(tdata,pRV,td_per);

pPAs  = smooth(td_per,pPAi, 20,'loess');
pPWs  = smooth(td_per,pPWi, 20,'loess');
pRAs  = smooth(td_per,pRAi, 20,'loess');
pRVs  = smooth(td_per,pRVi, 20,'loess');

% Shifting
pM_PA  = max(pPAs);
pM_RV  = max(pRVs);
pm_RA  = min(pRAs);
pm_PA  = min(pPAs);
dp     = pM_RV/1.01-pM_PA;
pPAs   = pPAs + dp;

shiftRV     = -5;  
pRV_noshift = pRVs;  
tempRV      = pRVs(end);
pRVs        = circshift(pRVs,shiftRV);
pRVs(end+shiftRV:end) = tempRV;

% Interpolate end point
midRV  = (pRVs(end+shiftRV-1)+pRVs(1))./[2];
midTRV = (td_per(end+shiftRV-1)+td_per(end))./[2];
pRVs(end+shiftRV-3:end) = interp1([td_per(end+shiftRV-3:end+shiftRV-1) midTRV td_per(end)],...
                                  [pRVs(end+shiftRV-3:end+shiftRV-1)' midRV  pRVs(1)],...
                                  td_per(end+shiftRV-3:end),'spline'); % Looks linear
                            
% Find when the max RV occurs
maxRV = find(pRVs == max(pRVs));
maxPA = find(pPAs == max(pPAs));

% Enforces that PA is one step ahead of RV
PAshift = (maxRV-maxPA)+1; 
pPAs = circshift(pPAs,PAshift);

% Shift the atrium in time
shiftRA = -8;    

% Determines the time from the unshifted signal, the shift those time values
[Tcrv, Trrv, Trra, tcra, Tcra] = AtrialTiming(pRAs, pRV_noshift, td_per); 
tempRA = pRAs(end);
pRAs   = circshift(pRAs,shiftRA);

% Timing parameters
Tcrv = 0.28;
Trrv = 0.51;

Trra = 0.15*T; 
tcra = 0.57; 
Tcra = 0.71; 

% Interpolate end point
try %If shift is negative this try will work 
midRA = (pRAs(end+shiftRA-1)+pRAs(1))./[2];
midTRA = (td_per(end+shiftRA-1)+td_per(end))./[2];

pRAs(end+shiftRA-3:end) = interp1([td_per(end+shiftRA-3:end+shiftRA-1) midTRA td_per(end)],...
                                  [pRAs(end+shiftRA-3:end+shiftRA-1)' midRA  pRAs(1)],...
                                   td_per(end+shiftRA-3:end),'spline'); % Looks linear
catch
end
% 
try
% Interpolate end point
%if shift is potive this try will work
midRA = (pRAs(end-shiftRA+1)+pRAs(1))./[2];
midTRA = (td_per(end-shiftRA+1)+td_per(end))./[2];

pRAs(end-shiftRA+3:end) = interp1([td_per(end-shiftRA+3:end-shiftRA+1) midTRA td_per(end)],...
                                  [pRAs(end-shiftRA+3:end-shiftRA+1)' midRA  pRAs(1)],...
                                   td_per(end-shiftRA+3:end),'spline'); % Looks linear
catch
end
                              
% Shift the atrium pressure up
dpRA = pRVs(1)-pRAs(1);
pRAs = pRAs + dpRA;
pRAs  = smooth(td_per,pRAs, 20,'loess');

% Try to set timing parameters from within patient file
data.tcra = tcra; 
data.Tcra = Tcra; 
data.Trrv = Trrv; 
data.Tcrv = Tcrv; 
data.Trra = Trra; 

% set up data file
data.T      = T;
data.NC     = NC;
data.dt     = dt;

data.td     = td;
data.td_per = td_per;

data.Gender = G;
data.Height = 161; 
data.Weight = 90; 
data.age    = age;
data.BSA    = 1.94;
data.CO     = 6.42;
data.pSAM   = 112; 
data.pSAm   = 76;
data.pSAa   = 1/3*data.pSAM+2/3*data.pSAm;
data.pPA    = pPAs;
data.pPW    = pPWs;
data.pRA    = pRAs;
data.pRV    = pRVs;

data.Patient = 1;

end

