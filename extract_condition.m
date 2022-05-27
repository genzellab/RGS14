%% CHANGE VARIABLE NAMES MANUALLY!
load('waveforms_cluster1_veh.mat');
data  = waveforms_cluster1_bp_veh;
HC_i = contains(data(:,8),'HC');
OS_i = ~contains(data(:,8),'HC');
waveforms_cluster1_bp_veh_HC = data(HC_i,:);
waveforms_cluster1_bp_veh_OS = data(OS_i,:);