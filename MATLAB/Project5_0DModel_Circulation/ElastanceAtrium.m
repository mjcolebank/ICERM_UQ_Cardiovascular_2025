function Ea = ElastanceAtrium(t,EaM,Eam,Tar,tac,Tac,T)
%Input the following variables:
%t = time from 0 to T
%T = length of cardiac cycle

if t<=Tar
   Ea = (EaM-Eam)*(1-cos(pi*(t-Tar)/(T-Tac+Tar)))/2+Eam;
elseif t <= tac    
  % Ea = (EaM-Eam)*sin(pi*(t-Tar)/(tac-Tar))/2+Eam;
   Ea = Eam; 
elseif t <= Tac
 %   Ea = (EaM-Eam)*sin(pi*(t-tac)/(Tac-tac))/2+Eam;
     Ea = (EaM-Eam)*(1-cos(pi*(t-tac)/(Tac-tac)))/2+Eam;
else
    Ea = (EaM-Eam)*(1+cos(pi*(t-Tac)/(T-Tac+Tar)))/2+Eam; 
end