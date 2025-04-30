function [IC, param] = get_CV_parameters()
% Now define the parameters for the system
P_LA     = 5; % left atrial pressure (mmHg)
P_SysCap = 20; % systemic capillary pressure (mmHg)
Rmv      = 5e-3; % resistance in the mitral valve (mmHg s/micro l)
Rav      = 1e-2; % resistance in the aortic valve (mmHg s/micro l)
Rart     = 1.2;  % resistance of the systemic arteries/arterioles (mmHg s/ml)
Cao      = 1.1;  % aortic compliance (ml s / mmHg)
Emax     = 1.3;  % End systolic elastance (mmHg/ml)
Emin     = 0.03;  % Nonlinear elastance term (1/ml) 
Vlv_d    = 5;    % volume not ejected in the heart
T_peak   = 0.25; % peak elastance (s)
T_relax  = 0.55;  % end of systole (s)
T        = 1.0;  % total cycle length

% Stack all the parameters into a vector
param = [P_LA,P_SysCap,Rmv,Rav,Rart,Cao,Emax,Emin,Vlv_d,T_peak,T_relax,T];

% Volumes are in microliters
Vlv_init = 200; % LV Volume
Vao_init = 380;  % Aortic Volume
IC = [Vlv_init; Vao_init]; % Our initial conditions
end