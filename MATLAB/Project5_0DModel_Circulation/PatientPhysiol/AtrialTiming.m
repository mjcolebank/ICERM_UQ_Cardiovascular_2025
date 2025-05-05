function [Tvcr,Tvrr,Tarr,tacr,Tacr] = AtrialTiming(pRA, pRV, tper)
% Timing parameters
% Note: assume we start from isovolumetric contraction of ventricles
% To find the timing parameters, we will take the time-series data, and
% determine when systole/diastole occur

% Find the average period length

%Find atrial relaxation parameters
diff_pRA = diff(pRA);

i=length(diff_pRA);
temp = -1;
while temp<0
    temp = diff_pRA(i);
    i=i-1;
end
Tacr  = tper(i);
pacrM = pRA(i);
    
while temp>0
    temp = diff_pRA(i);
    i=i-1;
end
tacr  = tper(i);
pacrm = pRA(i);


% Atrial contraction parameters
i = round(0.1*length(pRA));
temp = -1;
while temp<0
    temp = diff_pRA(i);
    i=i+1;
end
Tarr  = tper(i);
parr  = pRA(i);


% Ventricular timing parameters
diff_pRV = diff(pRV);
i = round(0.2*length(pRV));%20;
temp = 1;
while temp>=0
    temp = diff_pRV(i);
    i=i+1;
end
Tvcr = tper(i);
prvM = pRV(i);

i=i+10;
temp = -1;
while temp<=0
    temp = diff_pRV(i);
    i=i+1;
end
Tvrr = tper(i); 
prvm = pRV(i);


% Plot to ensure everything is kosher
% figure(2); clf;
% subplot(2,1,1);
% h = plot(tper,pRA,'r',Tacr,pacrM, '*k',tacr,pacrm,'og',Tarr,parr, 'sc');
% set(h(1),'Linewidth',3);
% set(h(2:4),'Markersize',10);
% 
% subplot(2,1,2); 
% h=plot(tper,pRV,'r',Tvrr,prvm,'*c',Tvcr,prvM,'*k');
% set(h(1),'Linewidth',3);
% set(h(2:3),'Markersize',10);

end