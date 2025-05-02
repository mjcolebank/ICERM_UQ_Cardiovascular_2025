function [pressure, flow] = ExtractData(nv,file_id)

% NOTE: Results are stored such that the first 1:N entries are the PROXIMAL
% large artery (1-15) and large vein (16-27) predictions, N+1:2N are the
% MIDPOINT predictions in the same vessels, and 2N+1:3N are the DISTAL
% predictions.

data = load(sprintf('output_%d.2d', file_id)); 

[~,~,p,q,~,~] = gnuplot(data);

pressure_all = (p(:,nv+1:2*nv)); % middle prediction
flow_all = (q(:,nv+1:2*nv)); % middle prediction

% real data consists of flow (entire time series) from vessels 1-7, 
% and pressure (min and max points only) from vessel 9, 
%%% so these are also the predictions we save
pressure  = pressure_all(:,9); %[max(pressure_all(:,9)), min(pressure_all(:,9))]; %???SHOULD THIS BE (:,9) OR (9,:)?
flow = flow_all(:,1:7); % vessel 8 is skipped

end