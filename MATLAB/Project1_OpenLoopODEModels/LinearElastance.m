% FUNCTION LinearElastance.m
% This is an introductory function example for the course on cardiovascular
% physiology: theory and practice.
%
% Author: Mitchel J. Colebank
% Date:5/17/2022

function p = LinearElastance(V,t,pars)
%Input the following variables:
%t = time from 0 to T
%T = heart rate

Emax    = pars(1); % mmHg/mL
Emin    = pars(2); % mmHg/mL
Vd      = pars(3); % mL
T_peak  = pars(4); % s
T_relax = pars(5); % s
T       = pars(6); % s

n = length(t);
E = zeros(1,n);
t = mod(t,T);
for i=1:n

    if t(i)<=T_peak
        E(i) = (Emax-Emin)*(1-cos(pi*t(i)/T_peak))/2 + Emin;
    elseif t(i) <= T_relax
        E(i) = (Emax-Emin)*(1+cos(pi*(t(i)-T_peak)/(T_relax-T_peak)))/2 + Emin;
    else
        E(i) = Emin;
    end
end
% Pressure is elastance times volume
p = E.*(V-Vd);

end
