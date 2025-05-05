% Calculate the axial and radial components of the tension vector

function [Tx,Ty] = get_subtensions(T,xm,ym)
sA = 2.*xm.*ym./(xm.^2+ym.^2);
cA = (-xm.^2 + ym.^2)./(xm.^2+ym.^2);
Tx = T.*sA;
Ty = T.*cA;
end