%--------------------------------------------------------------------------
%Used to fix some parameters and let the others vary (INDMAP) before
%solving ODE
%--------------------------------------------------------------------------

function [rout,J,CVsave] = model_wrap(pars,data)

tpars = data.ALLPARS;
tpars(data.INDMAP') = pars;

[rout,J,CVsave] = CVmodel(tpars,data);
end

