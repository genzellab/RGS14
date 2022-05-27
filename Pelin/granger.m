%% Granger Analysis

clc
clear
addpath('/Users/pelinozsezer/Documents/Science/Radboud/RGSProject/code/fieldtrip');
load('GC_ripple_4clusters_median_wa.mat');

%% What would you like to analyse?
input1 = GC_cluster1_veh_median_wa; % e.g., GC_cluster1_veh_median_wa
input2 = GC_cluster2_veh_median_wa; % e.g., GC_cluster2_veh_median_wa

input1_Bp = GC_Bp_cluster1_veh_median_wa; % e.g., GC_cluster1_veh_median_wa
input2_Bp = GC_Bp_cluster2_veh_median_wa; % e.g., GC_cluster2_veh_median_wa

% for normalization, add other cluster
input3 = GC_cluster3_veh_median_wa; % e.g., GC_cluster3_veh_median_wa

freqrange = [0:0.5:20]; % [0:0.5:20] or [20:1:100] or [100:2:300]

%% Compute granger
[granger,granger2,add]=granger_automation(input1,input2,input3,freqrange);

%% PLOTTING

    %% PFC-->HPC,2s       
    i=1;
    j=2;
    
    tf_p=squeeze(granger.grangerspctrm(i,j,:,:));
    tf_p2=squeeze(granger2.grangerspctrm(i,j,:,:)); 
    tf_p3=squeeze(add.grangerspctrm(i,j,:,:)); 
    
    %% Normalize the colorbar
    zmin= min([min(tf_p, [],'all'), min(tf_p2, [],'all'), min(tf_p3, [],'all')],[],'all');
    zmax= max([max(tf_p, [], 'all'), max(tf_p2, [],'all'), max(tf_p3, [],'all')],[],'all');
    clim_all =[zmin zmax];
    
    % input1   
    imagesc(-1:0.01:1,granger.freq,tf_p,clim_all); 
    axis xy % flip vertically
    colorbar
    colormap(hot(256))
    
    xlim([-1 1])
    xlabel('time')
    ylabel('frequency')
    title('PFC to HPC - Cluster 1 (Veh)')
    
    saveas(gcf,'pfc2hpc_veh_c1_020_2s.jpg');
    saveas(gcf,'pfc2hpc_veh_c1_020_2s.pdf');
    close all
    
    % input2      
    imagesc(-1:0.01:1,granger2.freq, tf_p2,clim_all); 
    axis xy % flip vertically
    colorbar
    colormap(hot(256))
    
    xlim([-1 1])
    xlabel('time')
    ylabel('frequency')
    title('PFC to HPC - Cluster 2 (Veh)') % change accordingly
    
    saveas(gcf,'pfc2hpc_veh_c2_020_2s.jpg'); % change accordingly
    saveas(gcf,'pfc2hpc_veh_c2_020_2s.pdf');
    close all
    
    % contrast        
    imagesc(-1.1:0.01:1.1,granger.freq,tf_p-tf_p2); 
    axis xy % flip vertically
    colorbar
    colormap(colorbar_cluster12) % change accordingly
    xlim([-1 1])
    xlabel('time')
    ylabel('frequency')
    title('PFC to HPC - Contrast: Cluster1-2 (Veh)'); % change accordingly
    
    saveas(gcf,'pfc2hpc_contrast_vehicle_c12_020_2s.jpg'); % change accordingly
    saveas(gcf,'pfc2hpc_contrast_vehicle_c12_020_2s.pdf');
    close all
    
    %% Stats
    
    %% input1
    p = input1;
    q = input1_Bp;
    
    %% Iteration of GC trials Veh 
    iter = 30;
    m = 400;
    grangerspctrm_concat = zeros(2,2,length(freqrange)-1,length([-1.1:0.01:1.1]),iter);
    for i = 1:iter
        i
        randorder = randperm(length(q));
        temp_q = q(randorder);
        temp_p = p(randorder);
        q_ = temp_q(1:m);
        p_ = temp_p(1:m);
        
        % Compute time frequency gc
        fn = 1000;
        leng = length(p_);
        ro = 3000;
        tm = create_timecell(ro,leng);
        label = [{'PFC'}; {'HPC'}];
        
        Data.label = label;
        Data.time = tm;
        Data.trial = p_.';
        
        cfg.bsfilter = 'yes';
        cfg.bsfreq = [49 51];
        Data = ft_preprocessing(cfg,Data); 
        
        [granger_tf] = createauto_timefreq(Data,freqrange,[-1.1:0.01:1.1]);
        grangerspctrm_concat(:,:,:,:,i) = granger_tf.grangerspctrm;
    end
    
    granger_tf= grangerspctrm_concat;   
    
    
    %% input2
    p = input2;
    q = input2_Bp;
    
    % Iteration of GC trials 
    iter = 30;
    m = 400;
    grangerspctrm_concat = zeros(2,2,length(freqrange)-1,length([-1.1:0.01:1.1]),iter);
    for i = 1:iter
        i
        randorder = randperm(length(q));
        temp_q = q(randorder);
        temp_p = p(randorder);
        q_ = temp_q(1:m);
        p_ = temp_p(1:m);
        
        % Compute time frequency gc
        fn = 1000;
        leng = length(p_);
        ro = 3000;
        tm = create_timecell(ro,leng);
        label = [{'PFC'}; {'HPC'}];
        
        Data.label = label;
        Data.time = tm;
        Data.trial = p_.';
        
        cfg.bsfilter = 'yes';
        cfg.bsfreq = [49 51];
        Data = ft_preprocessing(cfg,Data); 
        
        
        [granger_tf2] = createauto_timefreq(Data,freqrange,[-1.1:0.01:1.1]);
        grangerspctrm_concat2(:,:,:,:,i) = granger_tf2.grangerspctrm;
    end
    
    granger_tf2= grangerspctrm_concat2;
    
    a=1; %pfc
    b=2; %hpc
    
    [zmap]=stats_high_granger(granger_tf,granger_tf2,a,b);
    J=imagesc(granger.time,granger.freq,zmap);
    axis xy % flip vertically
    colorbar()
    colormap('colorbar_cluster12') % change accordingly
    J=title('PFC to HPC - Stats for Contrast: Cluster1-2 (Veh)'); % change accordingly
    J.FontSize=12;
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    xlim([-1 1])
    
    saveas(gcf,'pfc2hpc_stats_contrast_vehicle_c12_020_2s.jpg'); % change accordingly
    saveas(gcf,'pfc2hpc_stats_contrast_vehicle_c12_020_2s.pdf'); % change accordingly
    close all
    
    %% HPC->PFC 2s
    i=2;
    j=1;
    
    tf_p=squeeze(granger.grangerspctrm(i,j,:,:));
    tf_p2=squeeze(granger2.grangerspctrm(i,j,:,:));
    tf_p3=squeeze(add.grangerspctrm(i,j,:,:));
   
    zmin= min([min(tf_p2, [],'all'), min(tf_p, [],'all'), min(tf_p3, [],'all')],[],'all');
    zmax= max([max(tf_p2, [], 'all'), max(tf_p2, [],'all'), max(tf_p3, [],'all')],[],'all');
    clim =[zmin zmax];
    
    % input1  
    imagesc(-1:0.01:1,granger.freq,tf_p,clim_all); 
    axis xy % flip vertically
    colorbar
    colormap(hot(256))
    
    xlim([-1 1])
    xlabel('time')
    ylabel('frequency')
    title('HPC to PFC - Cluster 1 (Veh)') % change accordingly
    
    saveas(gcf,'hpc2pfc_veh_c1_020_2s.jpg'); % change accordingly
    saveas(gcf,'hpc2pfc_veh_c1_020_2s.pdf');
    close all
    
    % input2      
    imagesc(-1:0.01:1,granger2.freq, tf_p2,clim_all); 
    axis xy % flip vertically
    colorbar
    colormap(hot(256))
    
    xlim([-1 1])
    xlabel('time')
    ylabel('frequency')
    title('HPC to PFC - Cluster 2 (Veh)') % change accordingly
    
    saveas(gcf,'hpc2pfc_veh_c2_020_2s.jpg'); % change accordingly
    saveas(gcf,'hpc2pfc_veh_c2_020_2s.pdf');
    close all
       
    % contrast        
    imagesc(-1.1:0.01:1.1, granger.freq, tf_p-tf_p2); 
    axis xy % flip vertically
    colorbar
    colormap(colorbar_cluster12) % change accordingly
    xlim([-1 1])
    xlabel('time')
    ylabel('frequency')
    title('HPC to PFC - Contrast: Cluster1-2 (Veh)'); % change accordingly
    
    saveas(gcf,'hpc2pfc_contrast_vehicle_c12_020_2s.jpg'); % change accordingly
    saveas(gcf,'hpc2pfc_contrast_vehicle_c12_020_2s.pdf');
    close all
    
    %%stats
    
    %% input1
    p = input1;
    q = input1_Bp;
    
    %% Iteration of GC trials Veh 
    iter = 30;
    m = 400;
    grangerspctrm_concat = zeros(2,2,length(freqrange)-1,length([-1.1:0.01:1.1]),iter);
    for i = 1:iter
        i
        randorder = randperm(length(q));
        temp_q = q(randorder);
        temp_p = p(randorder);
        q_ = temp_q(1:m);
        p_ = temp_p(1:m);
        
        fn = 1000;
        leng = length(p_);
        ro = 3000;
        tm = create_timecell(ro,leng);
        label = [{'PFC'}; {'HPC'}];
        
        Data.label = label;
        Data.time = tm;
        Data.trial = p_.';
        
        cfg.bsfilter = 'yes';
        cfg.bsfreq = [49 51];
        Data = ft_preprocessing(cfg,Data); 
        
        [granger_tf] = createauto_timefreq(Data,freqrange,[-1.1:0.01:1.1]);
        grangerspctrm_concat(:,:,:,:,i) = granger_tf.grangerspctrm;
        
    end
    
    granger_tf= grangerspctrm_concat;   
    
    
    %%  input2
    p = input2;
    q = input2_Bp;
    
    %% Iteration of GC trials 
    iter = 30;
    m = 400;
    grangerspctrm_concat = zeros(2,2,length(freqrange)-1,length([-1.1:0.01:1.1]),iter);
    for i = 1:iter
        i
        randorder = randperm(length(q));
        temp_q = q(randorder);
        temp_p = p(randorder);
        q_ = temp_q(1:m);
        p_ = temp_p(1:m);
        
        % Compute time frequency gc
        fn = 1000;
        leng = length(p_);
        ro = 3000;
        tm = create_timecell(ro,leng);
        label = [{'PFC'}; {'HPC'}];
        
        Data.label = label;
        Data.time = tm;
        Data.trial = p_.';
        
        cfg.bsfilter = 'yes';
        cfg.bsfreq = [49 51];
        Data = ft_preprocessing(cfg,Data); 
        
        [granger_tf2] = createauto_timefreq(Data,freqrange,[-1.1:0.01:1.1]);
        grangerspctrm_concat2(:,:,:,:,i) = granger_tf2.grangerspctrm;
        
    end
    
    granger_tf2= grangerspctrm_concat2;
    
    % 2nd input is relative to that.
    a=2; %hpc
    b=1; %pfc
    
    [zmap]=stats_high_granger(granger_tf,granger_tf2,a,b);
    J=imagesc(granger.time,granger.freq,zmap);
    axis xy % flip vertically
    colorbar()
    colormap('colorbar_cluster12') % change accordingly
    J=title('HPC to PFC - Stats for Contrast: Cluster1-2 (Veh)'); % change accordingly
    J.FontSize=12;
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    xlim([-1 1])
    
    saveas(gcf,'hpc2pfc_stats_contrast_vehicle_c12_020_2s.jpg'); % change accordingly
    saveas(gcf,'hpc2pfc_stats_contrast_vehicle_c12_020_2s.pdf');
    close all