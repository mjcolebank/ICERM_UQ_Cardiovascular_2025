global ALLPARS INDMAP
%%%
% This script runs optimizations of subset 6 over various data shifts. This
% is ultimately a way to circumvent estimating initial conditions as
% parameters

tshift = -.2:.01:.2;
INDMAP = [2 8 9 10 11]; %parameters to optimize

D0 = data_process(0,2); %unshifted data
[x0,Init,low,hi] = load_global(D0); % Function returning x0 (the log of the parameters), Init (the initial conditions), 
%all of the variables loaded by load_global are the same for a given data shift
ALLPARS = x0;

%parameters and bounds
    optx   = x0(INDMAP); 
    opthi  = hi(INDMAP);
    optlow = low(INDMAP);

for i = 1:length(tshift)
    clear x xopt histout data
    
    data = data_process(tshift(i),2); %shifted data
    data.Init = Init; %pass IC into struct for CVmodel
    
    %optimizer settings
    maxiter = 30;                   % max number of iterations        
    mode    = 2;                    % Performs Levenberg-Marquart optimization
    nu0     = 2.d-1;                % Regularization parameter
    [xopt, histout] = newlsq_v2(optx,'opt_wrap',1.d-3,maxiter,mode,nu0,...
               opthi,optlow,data);
           
    %merging optimization results with fixed parameters       
    x = ALLPARS;                    % only optimized parameters
    x(INDMAP) = xopt;               % all parameters, including optimized one
    
    save(['shift_',num2str(tshift(i)),'_2.mat'])
end

exit;