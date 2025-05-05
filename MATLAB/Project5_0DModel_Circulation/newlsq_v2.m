function [x,histout,costdata, jachist, xhist,rout, sc,CVsave] = newlsq_v2(x0,f,tol,maxit,mode,nu0,high,low,data)

%
% C. T. Kelley, May 9, 2007
%
% This code comes with no guarantee or warranty of any kind.
%
% function [x,histout,costdata] = newlsq(x0,f,tol,maxit,mode,nu0,high,low)
%
% General Levenberg-Marquardt code: 
%         Your choice of linesearch, pseudo-transient continuation,
%         or (coming soon) trust region methods.
%
% Input: x0 = initial iterate
%        f  = rout^T rout/2 = objective function,
%            the calling sequence for f should be
%            [fout,gout,jac,rout]=f(x) where fout=f(x) is a scalar
%              gout = jac^T rout = grad f(x) is a COLUMN vector
%              and jac = rout' = Jacobian of rout is an M x N matrix
%              and rout=r;
%
% The residual is last in the argument list so I can give least squares
% problems to the other codes without using a wrapper.
%
%        tol = termination criterion norm(projgrad) < tol
%        maxit = maximum iterations (optional) default = 100
%        mode = 0 for pseudo-transient continuation,  SER-A timestep control
%               .1 for pseudo-transient continuation, SER-B timestep control
%               .2 for pseudo-transient continuation, TTE timestep control
%
%               3, 3.1, 3.2 for PTC with line search
%
%               1 for Levenberg-Marquardt with line search
%               2 for Levenberg-Marquardt with trust region
%               default = 0
%        nu0   = Levenberg parameter scaling
%                (mode = 0) nu0=1/dt0 
%                (mode = 1) nu=nu0*max(1,sqrt(f(x)))
%                default = 10 for mode 0, .1 for mode 1 and 2
%                Set nu0=0 to get Gauss-Newton.
%         high, low = vectors of upper and lower bounds. 
%
% Output: x = solution
%         histout = iteration history   
%         Each row of histout is      
%         [norm(projgrad), f, number of stepsize cuts, iteration count] 
%         costdata = [num f, num grad, num hess] (for levmar, num hess=0)
%
% At this stage all iteration parameters are hardwired in the code.
%
%

global treg y0 fe sc xb
fe = f;
highin = high;
lowin  = low;
x0in   = x0;
y0     = x0in;
%
% Set debug=0 to turn off the iteration-by-iteration stats.
%
debug=1;
%
% treg = Tikhonov regularization parameter.
%
treg=0.d0;
%
% hscale=0 does not scale the variables
% hscale=1 scales the variables by abs(x0)
% hscale=2 scales by the bounds
%
hscale=2;
%
% svalcut: The computation for the step will discard all singular values
%          < svalcut*(largest singular value)
%
svalcut=1.d-12;
%
% maxarm = limit on number of stepsize reductions in the line search
%
% numax  = max 1/dt or max Levenberg parameter.
%
maxarm=10;
numax =1.d4;
%
% Do not change lmversion. 
%
lmversion=2;
%
disp([mode, lmversion, hscale])
nc=length(x0);
[x0, sc, xb, high, low]=lsq_scale(highin, lowin, hscale, x0in);
itc=1; xc=x0; 

[fc,gc,jac,rout,CVsave]=feval(@tregf,xc,data);
jold=zeros(size(jac));
pgc=xc-proj(xc-gc,high,low);
mvar=length(gc); nvar=length(xc);
imode=round(mode);
ptcmode=mode-imode;
ptc_on=1;
switch imode
   case 0 % PTC + control of nu
      armijoflag=0;
      nu=nu0; nu0=nu0/norm(pgc);
   case 1 % LM with linesearch
      ptc_on=0;
      armijoflag=1;
      nu=nu0*min(1,sqrt(fc));
   case 2 % LM with TR (Higham)
      ptc_on=0;
      nu=nu0;
      armijoflag=2;
   case 3 % PTC-Armijo
      nu=nu0; nu0=nu0/norm(pgc);
      armijoflag=1;
   otherwise
      nu=nu0; nu0=nu0/norm(pgc);
      armijoflag=0;
end
nuc=nu; num=nuc; 
%
% Initialize quasi-Newton updates
%
B=zeros(nvar,nvar);
y=zeros(nvar,1);
step=zeros(nvar,1);
%
numf=1; numg=1; numh=0;
ithist(1,1)= norm(pgc); 
ithist(1,2)= fc; 
ithist(1,4)= itc-1; 
ithist(1,3)= 0;

%SRP
ithist(1,5) = cond(jac);
if nargout >= 4
    jachist{1} = jac;
    xhist{1}   = xb+sc*xc;
    %xct = xb+sc*xc;
    %x0t = xb+sc*x0;
    %xct
    %x0t
    %pause;
end

%SRPEND
npc=norm(pgc);
if debug == 1
   disp(ithist(itc,:))
end
armstop=0;
ared =1+tol*tol;
dcold=zeros(size(xc));
while(npc > tol & itc <= maxit & armstop == 0 & fc > tol*tol ...
       & abs(ared) > tol*tol)
        pact=proj_active(high,low,xc,npc);
        dc=moore(nu,jac,rout,svalcut,lmversion,pact);
%       [dc,B]=lsq_upd(nu,jac,rout,svalcut,lmversion,pact,jold,step,B);
        xt=proj(xc+dc,high,low); 
        dx=xt-xc;
        ddx=jac*dx;
%
%       Derive the correct formula for pred.
%
        pred=-(jac'*rout)'*dx-.5*ddx'*ddx;
%
        ft=feval(@tregf,xt,data);
        ared=fc-ft;
        if abs(ared) > tol*tol
        [xp, nup, idid, iflag]=lsq_new_point(xc,xt,nu,fc,gc,dc,ft,high,low,...
                               maxarm,ared,pred,armijoflag);
        else
           iflag=1;
           idid=0;
           nup=nu;
        end
        if iflag==-1
           armstop=1;
        end
        nu=nup;
        numf=numf+idid;
        if iflag==0
        stepm=step;
        step=xp-xc;
        xc=xp;
        gm=gc;
        jold=jac;
        fold=fc;
    	[fc,gc,jac,rout,CVsave]=feval(@tregf,xc,data); numf = numf+1; numg=numg+1;
        ared=fold-fc;
        y=gc-gm;
        pgc=xc-proj(xc-gc,high,low);
        npc=norm(pgc);
        itc=itc+1;
        ithist(itc,1)=norm(pgc); ithist(itc,2) = fc; 
        ithist(itc,4)=itc-1; ithist(itc,3)=idid;
        
        %SRP
        ithist(itc,5) = cond(jac);
        if nargout >= 4
            jachist{itc} = jac;
            xhist{itc}   = xb+sc*xc;
            %xct = xb+sc*xc;
            %x0t = xb+sc*x0;
            %xct
            %x0t
            %pause;
        end
        
        %SRPEND
        if debug == 1
             disp(ithist(itc,:))
        end
       end
        if ared > 0 | armijoflag == -1
        num=nuc;
        nuc=nu;
        if ptc_on == 1
           if itc == 1
           ptcmodex=ptcmode*(ptcmode~=3);
           nu=nuupdate(nu,step,nu0,npc,fc,ptcmodex,xc,num,stepm);
           else
           nu=nuupdate(nu,step,nu0,npc,fc,ptcmode,xc,num,stepm);
           end
        else
           nu=nuupdate(nu,step,nu0,npc,fc,mode);
        end
        else
        nu0=2.d0*nu0;
        end
% RENEE CAN YOU FIGURE THIS OUT        
%         if nu > numax 
%            disp('max nu exceeded');
%            armstop=1;
%         end
end
%x=xc; histout=ithist(1:itc,:);
x=xc; histout=ithist;
x=xb+sc*x;
costdata=[numf, numg, numh];


function [xp, nup, idid, iflag]=lsq_new_point(xc,xt,nu,fc,gc,dc,ft,high,low,...
                                              maxarm,ared,pred,lmode)
nup=nu; xp=xt; idid=0; iflag=0;
switch lmode;
case -1 % do nothing
case 0 % no line search; control nu based on reduction in f; 
        if ared < 0
           iflag=1;
           nup=2.d0*nu;
           xp=xc;
        end
case 1 % line search 
        [xp,idid,iflag]=lsq_armijo(xc,xt,fc,gc,dc,ft,high,low,maxarm);
case 2 % LM-TR
        [xp, nup, idid, iflag]=trust_regionb(xc, xt, nu, ared, pred);
otherwise
        disp('No new point rule given. Do nothing.');
end


function [xp,ididt,iflag]=lsq_armijo(xc,xt,fc,gc,dc,ft,high,low,maxarm)
        iflag=0;
        ididt=0; goal=fc-1.d-4*(xc-xt)'*(xc-xt);

        if ft <= goal 
                xp=xt;
        else
                [xp,idid]=polyline(xc,fc,gc,dc,ft,@tregf,maxarm,high,low);
                if idid==-1
                  disp('line search failure in newlsq.m');
                  xp=xc;
                  iflag=-1;
                end
                ididt=ididt+idid;
        end


function [xp, idid, lambda]=polyline(xc, fc, gc, d, ft, f, maxarm, high, low)
%
% C. T. Kelley, Dec 29, 1997
%
% This code comes with no guarantee or warranty of any kind.
%
% function [xp, idid]=polyline(xc, fc, gc, d, ft, fobj, maxarm, high, low)
%
% polynomial line search, call after first point is rejected
%
% Input: xc = current point
%        fc = current function value
%        gc = current gradient value
%         d = direction
%        ft = trial function (rejected value)
%         f = objective function
%             the calling sequence for f should be
%             [fout,gout]=f(x) where fout=f(x) is a scalar
%             and gout = grad f(x) is a COLUMN vector
%    maxarm = maximum number of step length reductions   
%    high,low = bounds on variables
%
% Output: xp = successful new point (if it exists)
%       idid = number of calls to f (if line search succeeds) or
%              -1 if line search fails.
%
% Requires: polymod.m
%
% line search parameters that everyone uses
%
alp=1.d-4; blow=.1; bhigh=.5;
%
% Set up the search
%
q0=fc; qp0=gc'*d; qc=ft; lamc=1; iarm=0; numf=0;
xt=proj(xc+d,high,low);
fgoal=fc-alp*lamc*(xc-xt)'*(xc-xt);
while ft > fgoal
    iarm=iarm+1;
    if iarm==1  % quadratic
       lambda=polymod(q0, qp0, lamc, qc, blow, bhigh);
    else
       lambda=polymod(q0, qp0, lamc, qc, blow, bhigh, lamm, qm);
    end
    qm=qc; lamm=lamc; lamc=lambda;
    xt=proj(xc+lambda*d,high,low);
    ft=feval(f,xt,data); 
    numf = numf+1; qc=ft;
    if(iarm > maxarm)
         disp(' line search failure'); idid=-1; xp=xc;
    return; end
    fgoal=fc-alp*lamc*(xc-xt)'*(xc-xt);
end
xp=xt; idid=numf;


function [lplus]=polymod(q0, qp0, lamc, qc, blow, bhigh, lamm, qm)
%
% C. T. Kelley, Dec 29, 1997
%
% This code comes with no guarantee or warranty of any kind.
%
% function [lambda]=polymod(q0, qp0, qc, blow, bhigh, qm)
%
% Cubic/quadratic polynomial linesearch
%
% Finds minimizer lambda of the cubic polynomial q on the interval
% [blow * lamc, bhigh * lamc] such that
%
% q(0) = q0, q'(0) = qp0, q(lamc) = qc, q(lamm) = qm
% 
% if data for a cubic is not available (first stepsize reduction) then
% q is the quadratic such that
% 
% q(0) = q0, q'(0) = qp0, q(lamc) = qc
%
lleft=lamc*blow; lright=lamc*bhigh; 
if nargin == 6
%
% quadratic model (temp hedge in case lamc is not 1)
%
    lplus = - qp0/(2 * lamc*(qc - q0 - qp0) );
    if lplus < lleft lplus = lleft; end
    if lplus > lright lplus = lright; end
else
%
% cubic model
%
    a=[lamc^2, lamc^3; lamm^2, lamm^3];
    b=[qc; qm]-[q0 + qp0*lamc; q0 + qp0*lamm];
    c=a\b;
    lplus=(-c(1)+sqrt(c(1)*c(1) - 3 *c(2) *qp0))/(3*c(2));
    if lplus < lleft lplus = lleft; end
    if lplus > lright lplus = lright; end
end


% Compute the search direction with an SVD. 
% Discard singular values < sfloor * (largest singular value).
% This function needs to know about the active set, which you
% tell it about via pact.
function dm=moore(nu,jac,rout,sfloor,lmversion,pact)
[mj,nj]=size(jac);
gc=jac'*rout;
pinact=eye(nj)-pact;
jac=jac*pinact;
[uj,sj,vj]=svd(jac,'econ');
if lmversion == 1 | nu < 1.d-12
   sd=diag(sj);
   sd=sd+1.d-14*sd(1);
   pj=(sd > sfloor*sd(1));
   nj=length(rout);
   sd=pj./sd;
   md=1+nu*sd.*sd;
   sd=sd./md;
%
% Version 2
% 
elseif lmversion == 2
   sd=diag(sj); pj=(sd > sfloor*sd(1)); 
   sd=pj.*sd;
   sd=sd./(nu+sd.*sd);
else
end
npro=1.d0/(1.d0+nu);
dm=-vj*diag(sd)*uj'*rout;
dm=pinact*dm-npro*pact*gc;


% Scaled and Tikhonov regularized function. The regularization parameter
% is a global within newlsq.m, so the calling function does not
% have to know about it.
%
% The scaling gets done by tregf, not in the main part of the code. 
function [fc,gc,jac,rout,CVsave]=tregf(zc,data)
global treg y0 fe sc xb
nc=length(y0);
xc=xb+sc*zc;
if nargout > 1
[fc,gc,jac,rout,CVsave]=feval(fe,xc,data);
jac=jac*sc;
gc=sc*gc;
if treg > 0
   fc=fc+.5*treg*(xc-y0)'*(xc-y0);
   rout=[rout', sqrt(treg)*(xc-y0)']';
   jac=[jac', sqrt(treg)*sc]';
   gc=jac'*rout;
end
else
fc=feval(fe,xc,data);
if treg > 0
   fc=fc+.5*treg*(xc-y0)'*(xc-y0);
end
end

% Projection onto the bounds.
function py=proj(y,high,low)
py=max(y,low); py=min(py,high);


% Compute the projection onto the epsilon-active indices.
% epsilon is computed in the standard way.
function pact=proj_active(high,low,xc,npc);
epsact =min(1.d-2,sqrt(npc));
limhigh=high-epsact*abs(high)-1.d-6*(high==0);
limlow =low+epsact*abs(low)+1.d-6*(low==0);
dact   =(xc < limlow) | (xc > limhigh);
pact   =diag(dact);

% Compute the scaling matrix and base vectors
function [x0, sc, xb, high, low]=lsq_scale(highin, lowin, hscale, x0in)
nc  = length(highin);
high= highin;
low = lowin;
x0  = x0in;
switch hscale
case 0
     sc=eye(nc);
     xb=zeros(nc,1);
case 1
     xb=zeros(nc,1);
     xl=(x0 ~=0); ml=min(abs(x0(xl))); dsc=(x0==0).*ml+abs(x0); 
     sc=diag(dsc);
case 2
     xb=lowin;
     sc=diag(highin-lowin);
otherwise
     error('hscale must be 0, 1, 2');
end
x0   = inv(sc)*(x0in-xb);     
high = inv(sc)*(highin-xb);
low  = inv(sc)*(lowin-xb);

% Update nu (or dti=1/dt)
function nu_out=nuupdate(nu_in,step,nu0,norm_projgrad,fval,umode,...
                         xc,num,stepm)
dumode=round(10*umode);
switch dumode
case 0 % PTC with SER-A
     nu_out=nu0*norm_projgrad;
case 1 % PTC with SER-B
     nu_out=nu_in*norm(step);
case 2 % PTC with TTE
     dtm=1/num; dtc=1/nu_in; 
     upp=(nu_in*step-num*stepm)*2/(dtc+dtm);
     upp=upp./(1+abs(xc));
     nu_out=sqrt(2*norm(upp,inf)/3);
case 10
     nu_out=nu0*min(1,sqrt(fval));
case 20
     nu_out=nu_in;
otherwise
     error(' must be 0, .1, .2, 2, 3, 3.1, 3.2');
end


function [xp, nup, idid, iflag]=trust_regionb(xc, xt, nu, ared, pred)
idid=0; iflag=0; xp=xt; nup=nu;
mulow=.25; muhigh=.75;
rat=ared/pred;
if ared < 0 
   xp=xc;
   nup=nu*2;
   iflag=1;
elseif rat < mulow
   xp=xt;
   nup=nu*2;
elseif rat <= muhigh & rat >=mulow
   xp=xt;
   nup=nu;
else
   xp=xt;
   nup=.5*nu;
end
