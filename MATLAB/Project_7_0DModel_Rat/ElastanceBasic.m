function Elv = ElastanceBasic(t,T,Tsf,Trf,Emin,EMax)
%Input the following variables:
%t = time from 0 to T
%T = heart rate
%Tsf = fraction of time to go from min to max elasticity
%Trf = fraction of time to go from max to min elasticity
%Em  = minimum elasticity
%EMm = maximum - minimum elasticity
%Outputs Elh (elasticity of the left ventricle)

Tr  = Trf*T;    % time it takes to go from max to min elasticity
Ts  = Tsf*T;  % time is takes to go from min to max elasticity

if t<=Ts
   Elv = Emin + ((EMax-Emin)/2)*(1-cos(pi*t/Ts)); 
elseif t <=Tr     
   Elv = Emin + ((EMax-Emin)/2)*(cos(pi*(t-Ts)/(Tr-Ts)) + 1);     
else
   Elv = Emin;
end