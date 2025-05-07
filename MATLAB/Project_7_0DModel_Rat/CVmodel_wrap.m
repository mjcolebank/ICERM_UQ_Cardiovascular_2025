function [J,rout] = CVmodel_wrap(pars,data)
% merges parameter for optimization w/ all model parameters
global ALLPARS INDMAP 

tpars = ALLPARS;
tpars(INDMAP') = pars;

[rout, J] = CVmodel(tpars,data); % for optimization
