% Get parameters postoptimization

clear; clc; close all;
Names = {'Rs','Rp','Rava','Rmva','Rpva','Rtva','Rpv','Rsv',... % 1-8
    'Csa','Csv','Cpa','Cpv',...                           % 9-12
    'EMra','Emra','EMla','Emla',...                       % 13-16
    'EMRv','EmRv','EMLv','EmLv',...                       % 17-20
    'Trra','tcra','Tcra','Tcrv','Trrv'};                  % 21-25

all_pars = zeros(9,8,25);
R=2;
for ptnb = 1:9
    for iter=1:8
        s = strcat('Opt_',num2str(ptnb),'_R',num2str(R),'_iter',num2str(iter),'.mat');

        cd OUT
        load(s);
        pars = exp(optpars);
        cd ..
        all_pars(ptnb,iter,:) = pars;
    end
end

%%
mean_pars = squeeze(mean(all_pars,2));
sd_pars   = squeeze(std(all_pars,[],2));

CoV = sd_pars(:,data.INDMAP)./mean_pars(:,data.INDMAP);

h=figure;
h=boxplot(CoV)
set(h,'linewidth',2,'markersize',10)
set(gca,'FontSize',18);
title('Res 1: All pars, noTrra, noEmla, noCpv');
xticklabels(Names(data.INDMAP))
xtickangle(45);
grid on;
ylabel('Coeff Var')
ax = gcf;
exportgraphics(ax,'Boxplot.png','Resolution',300)
