R = 2; % 1 static only & 2 static & dynamic

ind = {};
for ptnb = 1
     INDMAP  = [1 2 6 8 9:11 14    18 20 22:25]; 
     IND     = [1 2 6 8 9:11 14    18 20]; 
     
     % x0 = [Rs Rp Rava Rmva Rpva Rtva Rpv Rsv ...  % 1-8
     %      Csa Csv Cpa Cpv ...                     % 9-12
     %      EMra Emra EMla Emla ...                 % 13-16
     %      EMrv Emrv EMlv Emlv ...                 % 17-20
     %      Trra tcra Tcra ...                      % 21-23
     %      Tcrv Trrv]';                            % 24-25
          
     %task = 1;     % Run Sens
     %task = 2;     % Run Opt
      task = 3;     % Plot Results (Run's forward model)
     %task = 4;     % Plot Sens
        
     par = 1; % Nom pars org (0), Opt pars (2)
     display(strcat('patient ',num2str(ptnb)));
     rng('Shuffle');
     rnum  = zeros(25,9);
     rt = -1+2*rand(length(IND),8);
     rnum(IND,1:8) = rt;
     iternum = [1 2 3 4 5 6 7 8 9];
     for i = 1:1
         i
         iter = iternum(i);
         RunSimulation(R,INDMAP,task,par,ptnb,iter,rnum(:,iter));  
     end
    disp('end patient');
    %pause;
    clear IND INDMAP i iter iternum par rnum rt task
    %close all
end

