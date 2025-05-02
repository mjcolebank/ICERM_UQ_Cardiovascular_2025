function [pressure, flow, pass] = Run_fluidsModel(x, file_id, nv, ntp)

f2_vec = x(:,1);
f3_vec = x(:,2);
fs2_vec = x(:,3);
fs3_vec = x(:,4);
alpha_vec = x(:,5);

n = size(x,1);

pass = NaN(n,1);

for i = 1:n
    
    % call the model
    param = [f2_vec(i), f3_vec(i), fs2_vec(i), fs3_vec(i), alpha_vec(i), file_id];
    
    param_str = mat2str(param);
    
    % Run the model
    % NOTE: Windows users need 'sor06.exe', Mac/Linux users need ./sor06
    cx = unix(sprintf('./sor06  %s',param_str(2:end-1)));
    
    if cx == 1 % marks a successful simulation
        
        pass(i) = 1;
        
        [pressure, flow] = ExtractData(nv,file_id);
        
    else
        
        pass(i) = 0;
        
        flow = NaN(ntp,nv-2);
        pressure = NaN(ntp,1);
    end
end