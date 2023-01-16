
%% SPECTROGRAM
clc
clear
addpath('/Users/pelinozsezer/Documents/Science/Radboud/RGS-Project/code/ripples/fieldtrip');
addpath('/Users/pelinozsezer/Documents/Science/Radboud/RGS-Project/code/ripples/functions');
load('GC_ripple_4clusters_median_wa.mat');

%% Merge data per treatment
GC_veh=[GC_Bp_cluster1_veh_median_wa; GC_Bp_cluster2_veh_median_wa;...
    GC_Bp_cluster3_veh_median_wa; GC_Bp_cluster4_veh_median_wa];

GC_rgs=[GC_Bp_cluster1_rgs_median_wa; GC_Bp_cluster2_rgs_median_wa;...
    GC_Bp_cluster3_rgs_median_wa; GC_Bp_cluster4_rgs_median_wa];

%% What would you like to analyse?
input1 = GC_veh;
input2 = GC_rgs;

channel   = 'PFC'; % 'HPC' or 'PFC'
freqrange = [0:0.5:20]; % [0:0.5:20] or [20:1:100] or [100:2:300]

%% Compute spectrogram
[freq,freq2,subs_freq,zlim] = spectrogram_automation(input1,input2,channel,freqrange);

%% PLOTTING

% input1
cfg              = [];
cfg.zlim         = zlim;
cfg.channel      = channel;
cfg.colormap     = colormap(hot);

ft_singleplotTFR(cfg, freq); 
g=title('PFC - Veh'); 
g.FontSize=12;
xlabel('Time (s)')
ylabel('Frequency (Hz)')
xlim([-1 1])

saveas(gcf,'020_2s_veh.jpg'); 
saveas(gcf,'020_2s_veh.pdf');
close all

% input2
cfg              = [];
cfg.zlim         = zlim;
cfg.channel      = channel;
cfg.colormap     = colormap(hot);

ft_singleplotTFR(cfg, freq2); 
g=title('PFC - RGS'); 
g.FontSize=12;
xlabel('Time (s)')
ylabel('Frequency (Hz)')
xlim([-1 1])

saveas(gcf,'020_2s_rgs.jpg');
saveas(gcf,'020_2s_rgs.pdf');
close all

% Contrast
cfg              = [];
cfg.channel      = channel;

ft_singleplotTFR(cfg, subs_freq); 
colormap(colorbar_treatment);
g=title('PFC - Contrast: Veh-RGS'); 
g.FontSize=12;
xlabel('Time (s)')
ylabel('Frequency (Hz)')
xlim([-1 1])

saveas(gcf,'contrast_veh-rgs_020_2s.jpg');
saveas(gcf,'contrast_veh-rgs_020_2s.pdf');


%% STATS

% Relative to the second input.
zmap=stats_high_spec(freq2,freq,1); % PFC=1 & HPC=2; output= freq x time

zmap(zmap == 0) = NaN; % convert 0s to nans in zmap - 0s are significantly not different
J=imagesc(freq.time,freq2.freq,zmap)
xlabel('Time (s)'), ylabel('Frequency (Hz)')
set(gca,'xlim',xlim,'ydir','no')
set(J,'AlphaData',~isnan(zmap))
colorbar()
colormap('colorbar_treatment') 

J=title('PFC - Stats for Contrast: Veh-RGS');
J.FontSize=12;
xlabel('Time (s)')
ylabel('Frequency (Hz)')
xlim([-1 1])

saveas(gcf,'stats_contrast_veh-rgs_020_2s.jpg'); 
saveas(gcf,'stats_contrast_veh-rgs_020_2s.pdf');
close all