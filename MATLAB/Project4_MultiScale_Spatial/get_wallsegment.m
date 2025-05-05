% This function takes the current midwall radius (y) and cavity volume and
% produce the distance from the wall to the cavity
% Based on the numerical solution from mathematica
function [xm, Am, Cm] = get_wallsegment(Vm,ym)
Vtemp = Vm.*(3./pi);
% Q = (Vtemp.^2 + sqrt(Vtemp.^2+ym.^6)).^(1/3); %old
Q = (Vtemp + sqrt(Vtemp.^2+ym.^6)).^(1/3); %new

xm = Q-ym.^2./Q;

Am = pi.*(xm.^2+ym.^2);
Cm = 2.*xm./(xm.^2+ym.^2);

end