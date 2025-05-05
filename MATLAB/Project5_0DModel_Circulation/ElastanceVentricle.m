function Elv = ElastanceVentricle(t,EvM,Evm,Tvc,Tvr)
%Input the following variables:
%t = time from 0 to T
%T = heart rate

if t<=Tvc
   Elv = (EvM-Evm)*(1-cos(pi*t/Tvc))/2 + Evm;
elseif t <= Tvr
   Elv = (EvM-Evm)*(1+cos(pi*(t-Tvc)/(Tvr-Tvc)))/2 + Evm;
else 
   Elv = Evm;
end

