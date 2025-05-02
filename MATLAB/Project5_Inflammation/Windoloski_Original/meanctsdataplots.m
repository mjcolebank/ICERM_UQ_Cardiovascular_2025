function meanctsdataplots(i,plotpars)

    global dataC

    load 'CTSmeanSDNOOutlierdata.mat';
    
    dataC.TNF  = TNFm;
    dataC.IL6  = IL6m;
    dataC.IL8  = IL8m;
    dataC.IL10 = IL10m;
    dataC.time = [0 1 2 3 4 6];
    
    height = plotpars.height;
    width  = plotpars.width;
    lw     = plotpars.lw;
    time = [0 1 2 3 4 6];

    figure(i);
    subplot(2,4,5); hold on;
    errorbar(time,TNFm,TNFsd,'ro','linewidth',lw);
    ylabel('TNF-a (pg/mL)');
    legend('');
    fig=gcf;
    fig.Position(3:4)=[width,height];
    
    subplot(2,4,6); hold on;
    errorbar(time,IL6m,IL6sd,'ro','linewidth',lw);
    ylabel('IL-6 (pg/mL)');
    legend('');
    fig=gcf;
    fig.Position(3:4)=[width,height];
    
    subplot(2,4,7); hold on;
    errorbar(time,IL8m,IL8sd,'ro','linewidth',lw);
    xlabel ('Time (hours)');
    ylabel('IL-8 (pg/mL)');
    legend('');
    fig=gcf;
    fig.Position(3:4)=[width,height];
    
    subplot(2,4,8); hold on;
    errorbar(time,IL10m,IL10sd,'ro','linewidth',lw);
    xlabel('Time (hours)');
    ylabel('IL-10 (pg/mL)');
    legend('');
    fig=gcf;
    fig.Position(3:4)=[width,height];

end