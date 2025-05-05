function [history, x] = DriverBasic_opt(data,pars)
% x0 = [Rs Rp Rava Rmva Rpva Rtva Rpv Rsv ...  % 1-8
%       Csa Csv Cpa Cpv ...                    % 9-12
%       EMra Emra EMla Emla ...                % 13-16
%       EMrv Emrv EMlv Emlv ...                % 17-20
%       Trra tcra Tcra Tcrv Trrv]';            % 21-25

%par_cluster = 0;

gPars.ALLPARS = pars;   % all parameters
gPars.INDMAP  = data.INDMAP; 

optx   = pars(gPars.INDMAP);

%set up structure for optimization
data.ALLPARS  = pars;

opthi     = data.hi(data.INDMAP);
optlow    = data.low(data.INDMAP);

% LM Kelley
% optimization information
maxiter = 50;        % max number of iterations        
mode    = 2;         % Performs Levenberg-Marquart optimization
nu0     = 2.d-1;     % Regularization parameter
[xopt, histout, costdata, jachist, xhist, rout, sc] = ...
    newlsq_v2(optx,'opt_wrap',1.d-3,maxiter,mode,nu0,...
              opthi,optlow,data);
history.x      = xhist;
history.hist   = histout;

% max_fun_evals = 1000;
% history = struct('fval',[],'x',[],'gradient',[]);
% 
% options = optimoptions(@fmincon,'Display','iter',...
%        'OutputFcn',@outfun,'MaxIter',maxiter,'MaxFunEvals',max_fun_evals,...
%        'FunctionTolerance',1e-3,'StepTolerance',0.001,'Display','iter');   
% 
% problem = createOptimProblem('fmincon','objective',...
%          @(optx)opt_wrap(optx,data),...
%          'x0',optx,'lb',optlow,'ub',opthi,'options',options);
%      
% [xopt,fval,exitflag,output,lambda,grad,hessian] = fmincon(problem);
% 
% 
% history.xopt     = xopt;
% history.exitflag = exitflag;
% history.output   = output;
% history.lambda   = lambda;
% history.grad     = grad;
% history.hessian  = hessian;
% history.fval     = fval;

history.INDMAP   = data.INDMAP;

x = data.ALLPARS;
x(data.INDMAP) = xopt;
   
% function stop = outfun(x,optimValues,state)
%      stop = false;
%  
%      switch state
%          case 'init'
%              hold on
%          case 'iter'
%          % Concatenate current point and objective function value with history
%          % x must be a row vector
%          history.fval = [history.fval; optimValues.fval];
%          history.x = [history.x; x];
%          history.gradient = [history.gradient; optimValues.gradient];
%           
%          case 'done'
%              hold off
%          otherwise
%      end
%end
end


