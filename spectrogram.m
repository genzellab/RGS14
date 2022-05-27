%% SPECTROGRAM
clc
clear
addpath('/Users/pelinozsezer/Documents/Science/Radboud/RGSProject/code/fieldtrip');
load('GC_ripple_4clusters_median_wa.mat');

%% What would you like to analyse?
input1 = GC_cluster1_veh_median_wa; % e.g., GC_cluster1_veh_median_wa
input2 = GC_cluster2_veh_median_wa; % e.g., GC_cluster2_veh_median_wa

% for normalization, add other cluster
input3 = GC_cluster3_veh_median_wa; % e.g., GC_cluster3_veh_median_wa

channel   = 'HPC'; % 'HPC' or 'PFC'
freqrange = [0:0.5:20]; % [0:0.5:20] or [20:1:100] or [100:2:300]

%% Compute spectrogram
[freq,freq2,subs_freq,zlim] = spectrogram_automation(input1,input2,input3,channel,freqrange);

%% PLOTTING

% input1
cfg              = [];
cfg.zlim         = zlim;
cfg.channel      = channel;
cfg.colormap     = colormap(hot);

ft_singleplotTFR(cfg, freq); 
g=title('HPC - Veh: Cluster 1'); % change accordingly
g.FontSize=12;
xlabel('Time (s)')
ylabel('Frequency (Hz)')
xlim([-1 1])

saveas(gcf,'cluster1_020_2s_veh.jpg'); % change accordingly
saveas(gcf,'cluster1_020_2s_veh.pdf');
close all

% input2
cfg              = [];
cfg.zlim         = zlim;
cfg.channel      = channel;
cfg.colormap     = colormap(hot);

ft_singleplotTFR(cfg, freq2); 
g=title('HPC - Veh: Cluster 2'); % change accordingly
g.FontSize=12;
xlabel('Time (s)')
ylabel('Frequency (Hz)')
xlim([-1 1])

saveas(gcf,'cluster2_020_2s_veh.jpg');
saveas(gcf,'cluster2_020_2s_veh.pdf');
close all


% Contrast
cfg              = [];
cfg.channel      = channel;

ft_singleplotTFR(cfg, subs_freq); 
colormap(colorbar_cluster12); % change accordingly
g=title('HPC - Contrast: Vehicle (Cluster1-2)'); % change accordingly
g.FontSize=12;
xlabel('Time (s)')
ylabel('Frequency (Hz)')
xlim([-1 1])

saveas(gcf,'contrast_vehicle_c12_020_2s.jpg'); % change accordingly
saveas(gcf,'contrast_vehicle_c12_020_2s.pdf');


%% STATS

% 2nd input is relative to that.
zmap=stats_high_spec(freq2,freq,2); % PFC=1 & HPC=2; output= freq x time

zmap(zmap == 0) = NaN; % convert 0s to nans in zmap - 0s are significantly not different
J=imagesc(freq.time,freq2.freq,zmap)
xlabel('Time (s)'), ylabel('Frequency (Hz)')
set(gca,'xlim',xlim,'ydir','no')
set(J,'AlphaData',~isnan(zmap))
colorbar()
colormap('colorbar_cluster12') % change accordingly

J=title('HPC - Stats for Contrast: Vehicle (Cluster1-2)');
J.FontSize=12;
xlabel('Time (s)')
ylabel('Frequency (Hz)')
xlim([-1 1])

saveas(gcf,'stats_contrast_vehicle_c12_020_2s.jpg'); % change accordingly
saveas(gcf,'stats_contrast_vehicle_c12_020_2s.pdf');
close all

