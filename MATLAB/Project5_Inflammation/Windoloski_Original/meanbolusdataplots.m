function meanbolusdataplots(i,plotpars)

    global dataB

    load 'BolusDataCorrectedSD.mat';
    
    TNFm = TNFm_NoOutliers;
    IL6m = IL6m_NoOutliers;
    IL8m = IL8m_NoOutliers;
    IL10m = IL10m_NoOutliers;
    
    dataB.TNF  = TNFm;
    dataB.IL6  = IL6m;
    dataB.IL8  = IL8m;
    dataB.IL10 = IL10m;
    dataB.time = Time;
    
    height = plotpars.height;
    width  = plotpars.width;
    lw     = plotpars.lw; 
    
    time = Time;

    figure(i);
    subplot(2,4,5); hold on;
    errorbar(time,TNFm,TNFsd,'k^','linewidth',lw);
    ylabel('TNF-a (pg/mL)');
    legend('');
    fig=gcf;
    fig.Position(3:4)=[width,height];
    
    subplot(2,4,6); hold on;
    errorbar(time,IL6m,IL6sd,'k^','linewidth',lw);
    ylabel('IL-6 (pg/mL)');
    legend('');
    fig=gcf;
    fig.Position(3:4)=[width,height];
    
    subplot(2,4,7); hold on;
    errorbar(time,IL8m,IL8sd,'k^','linewidth',lw);
    xlabel ('Time (hours)');
    ylabel('IL-8 (pg/mL)');
    legend('');
    fig=gcf;
    fig.Position(3:4)=[width,height];
    
    subplot(2,4,8); hold on;
    errorbar(time,IL10m,IL10sd,'k^','linewidth',lw);
    xlabel('Time (hours)');
    ylabel('IL-10 (pg/mL)');
    legend('');
    fig=gcf;
    fig.Position(3:4)=[width,height];

end