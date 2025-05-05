function [J, gout, jac, rout,CVsave] = opt_wrap(pars,data)
% note pars ONLY contains parameters to be estimated!
% function returns output for optimizer (J, gout, jac, rout)

[rout,J,CVsave]  = model_wrap(pars,data);
[jac]            = diffjac(pars,@myfun,rout,data);  
gout             = jac' * rout;


function rout = myfun(pars,data)

[rout,J] = model_wrap(pars,data);



function [jac] = diffjac(x,f,f0,data)
% Compute a forward difference dense Jacobian f'(x), return lu factors.
%
% uses dirder
%
% C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%

% inputs:
%          x, f = point and function
%          f0   = f(x), preevaluated

n = length(x);
jac = zeros(length(f0), n);
for j = 1:n
    zz = zeros(n,1);
    zz(j) = 1;
    jac(:,j) = dirder(x,zz,f,f0,data);
end


function z = dirder(x,w,f,f0,data)
% Compute a finite difference directional derivative.
% Approximate f'(x) w
% 
% C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
% function z = dirder(x,w,f,f0)
%
% inputs:
%           x, w = point and direction
%           f = function
%           f0 = f(x), in nonlinear iterations
%                f(x) has usually been computed
%                before the call to dirder

% Hardwired difference increment.
epsnew = data.DIFF_INC;

n = length(x);

% scale the step
if norm(w) == 0
    z = zeros(n,1);
return
end

% Now scale the difference increment.
xs=(x'*w)/norm(w);
if xs ~= 0.d0
    epsnew=epsnew*max(abs(xs),1.d0)*sign(xs);
end
epsnew=epsnew/norm(w);

% del and f1 could share the same space if storage
% is more important than clarity.
del = x+epsnew*w;
f1  = feval(f,del,data);
z   = (f1 - f0)/epsnew;


% Compute the step length with the three point parabolic model.
function [step,iarm,xp,fp,armflag] = armijo(direction,x,f0,f,maxarm)
iarm   = 0;
sigma1 = .5;
alpha  = 1.d-4;
armflag = 0;
xp = x; fp = f0; 

xold   = x;
lambda = 1; lamm = 1; lamc = lambda; iarm = 0;
step = lambda*direction;
xt   = x + step;
ft   = feval(f,xt);
nft  = norm(ft); nf0 = norm(f0); ff0 = nf0*nf0; ffc = nft*nft; ffm = nft*nft;
while nft >= (1 - alpha*lambda) * nf0;

    %   Apply the three point parabolic model.
    if iarm == 0
         lambda = sigma1*lambda;
    else
         lambda = parab3p(lamc, lamm, ff0, ffc, ffm);
    end

    % Update x; keep the books on lambda.
    step = lambda*direction;
    xt   = x + step;
    lamm = lamc;
    lamc = lambda;

    % Keep the books on the function norms.
    ft   = feval(f,xt);
    nft  = norm(ft);
    ffm  = ffc;
    ffc  = nft*nft;
    iarm = iarm+1;
    if iarm > maxarm
        disp(' Armijo failure, too many reductions ');
        armflag = 1;
        sol = xold;
        return;
    end
end
xp = xt; fp = ft;
% end of line search


function lambdap = parab3p(lambdac, lambdam, ff0, ffc, ffm)
% Apply three-point safeguarded parabolic model for a line search.
%
% C. T. Kelley, April 1, 2003
%
% This code comes with no guarantee or warranty of any kind.
%
% function lambdap = parab3p(lambdac, lambdam, ff0, ffc, ffm)
%
% input:
%       lambdac = current steplength
%       lambdam = previous steplength
%       ff0 = value of \| F(x_c) \|^2
%       ffc = value of \| F(x_c + \lambdac d) \|^2
%       ffm = value of \| F(x_c + \lambdam d) \|^2
%
% output:
%       lambdap = new value of lambda given parabolic model
%
% internal parameters:
%       sigma0 = .1, sigma1 = .5, safeguarding bounds for the linesearch

% Set internal parameters.
sigma0 = .1; sigma1 = .5;

% Compute coefficients of interpolation polynomial.
%
% p(lambda) = ff0 + (c1 lambda + c2 lambda^2)/d1
%
% d1 = (lambdac - lambdam)*lambdac*lambdam < 0
%      so, if c2 > 0 we have negative curvature and default to
%      lambdap = sigam1 * lambda.
c2 = lambdam*(ffc-ff0)-lambdac*(ffm-ff0);
if c2 >= 0
    lambdap = sigma1*lambdac; return
end
c1 = lambdac*lambdac*(ffm-ff0)-lambdam*lambdam*(ffc-ff0);
lambdap = -c1*.5/c2;
if lambdap < sigma0*lambdac, lambdap = sigma0*lambdac; end
if lambdap > sigma1*lambdac, lambdap = sigma1*lambdac; end
