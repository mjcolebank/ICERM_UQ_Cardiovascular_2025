function data = data_process(tshift,i)
tshift

%%%
% For the purposes of this work we are only interested in the basline time 
% series information prior to blood withdrawal. This function extracts a 
% baseline segments of pressure and volume data from rat blood withdrawal 
% experiments. This function shifts the time point at which the
% data beginnings in order to mitigate challenges with aligning the phase
% of data and model.
%%%%

s = strcat('DataBW',num2str(i),'.mat');
load(s); %Raw data - sampled at 1 kHz and contains blood withdrawal
% t - time vector (sec)
% P - pressure time series (mmHg)
% V - volume time series (uL)
% t_per - time stampts that denote the beginnning of systole
% PMax - maximal pressure values for each cardiac cycle
% Pmin - minimal pressure values for each cardiac cycle
% VMax - maximal volume values for each cardiac cycle
% Vmin - minimal volume values for each cardiac cycle

%%%
EPS = 1e-6;  
dt  = mean(diff(t)); 
N   = round(tshift/dt);
Ip1 = 2; 
tp1 = t_per(Ip1);
Ip2 = find(t_per<19+EPS, 1, 'last' );
tp2 = t_per(Ip2);
data.t_per = t_per(Ip1:Ip2)-t_per(Ip1);
data.V_per = V_per(Ip1:Ip2);
ID1 = find(t<tp1+eps, 1, 'last' );
ID2 = find(t<tp2+eps, 1, 'last' );
data.t = t(ID1:ID2)-t(ID1);
data.V = V(ID1+N:ID2+N);
data.P = P(ID1+N:ID2+N);
data.VMax = VMax(Ip1:Ip2);
data.Vmin = Vmin(Ip1:Ip2);
data.PMax = PMax(Ip1:Ip2);
data.Pmin = Pmin(Ip1:Ip2);


