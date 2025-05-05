%--------------------------------------------------------------------------
%Computes the sensitivity matrix dy/dpars.
%For relative sensitivity matrix dy/dlog(pars) = dy/dpars*pars, set 
%pars = log(pars) in DriverBasic_sens.
%--------------------------------------------------------------------------

function [sens] = senseq_cs(pars,data)

h = data.DIFF_INC;

for j = 1:length(pars)
    
    pars1    = pars;
    pars1(j) = pars(j) + sqrt(-1)*h;
    
    [y_im] = CVmodel(pars1,data);
    sens(:,j) = imag(y_im)./h;
end






    
