%--------------------------------------------------------------------------
%Computes the sensitivity matrix dy/dpars.
%For relative sensitivity matrix dy/dlog(pars) = dy/dpars*pars, set
%pars = log(pars) in DriverBasic_sens
%Plots ranked sensitivities
%--------------------------------------------------------------------------
function history = DriverBasic_sens(data,pars)

data.Names = {'Rs','Rp','Rava','Rmva','Rpva','Rtva','Rpv','Rsv',...
         'Csa','Csv','Cpa','Cpv',...
         'EMRa','EmRa','EMLa','EmLa',...
         'EMRv','EmRv','EMLv','EmLv',...
         'Trra','tcra','Tcra','Tcrv','Trrv'};

%senseq finds the non-weighted sensitivities 
[sens] = senseq_cs(pars,data);

% ranked classical sensitivities
[M,N] = size(sens);
for i = 1:N
    sens_norm(i)=norm(sens(:,i),2);
end

[Rsens,Isens] = sort(sens_norm,'descend');

history.Rsens = Rsens;
history.Isens = Isens;
history.sens  = sens;
history.Nsens = sens_norm;
history.data  = data;

end
