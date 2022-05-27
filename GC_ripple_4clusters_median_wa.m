%% MANUAL SCRIPT!!!
clc
clear
load(''); % all waveforms of the clusters

%% Prepare data for Granger Analysis!
% cluster 1
GC_cluster1_veh           = waveforms_cluster1_raw_veh(:,1);
GC_Bp_cluster1_veh        = waveforms_cluster1_bp_veh(:,1);
GC_time_cluster1_veh     = waveforms_cluster1_raw_veh(:,2:4);

% cluster 2
GC_cluster2_veh           = waveforms_cluster2_raw_veh(:,1);
GC_Bp_cluster2_veh        = waveforms_cluster2_bp_veh(:,1);
GC_time_cluster2_veh      = waveforms_cluster2_raw_veh(:,2:4);

% cluster 3
GC_cluster3_veh           = waveforms_cluster3_raw_veh(:,1);
GC_Bp_cluster3_veh        = waveforms_cluster3_bp_veh(:,1);
GC_time_cluster3_veh      = waveforms_cluster3_raw_veh(:,2:4);

% cluster 4
GC_cluster4_veh           = waveforms_cluster4_raw_veh(:,1);
GC_Bp_cluster4_veh        = waveforms_cluster4_bp_veh(:,1);
GC_time_cluster4_veh      = waveforms_cluster4_raw_veh(:,2:4);

%% VEH

%% per cluster 
R = (cellfun(@(equis1) max(abs(hilbert(equis1(2,3001-50:3001+50)))),GC_Bp_cluster1_veh));

[~,r_nl]=sort(abs(R-median(R)),'ascend');
R=R(r_nl);
GC_cluster1_veh=GC_cluster1_veh(r_nl);
GC_cluster1_veh = GC_cluster1_veh(1:2011);

GC_Bp_cluster1_veh=GC_Bp_cluster1_veh(r_nl);
GC_Bp_cluster1_veh= GC_Bp_cluster1_veh(1:2011);

GC_time_cluster1_veh=GC_time_cluster1_veh(r_nl,:);
GC_time_cluster1_veh = GC_time_cluster1_veh(1:2011,:);

fn=1000;
leng=length(GC_cluster1_veh);
ro=3000;
tm = create_timecell(ro,leng);

label=[{'PFC'}; {'HPC'}];
Data.label=label;
Data.time=tm;
Data.trial=GC_cluster1_veh.';

% NOTCH
cfg.bsfilter = 'yes';
cfg.bsfreq = [49 51];
Data = ft_preprocessing(cfg,Data);

% Run bandpass filter for plotting
cfg.bpfilter = 'yes';
cfg.bpfreq = [100 300];
Data = ft_preprocessing(cfg,Data);

% Detect artifact (overlay all ripples)
x=[];
for i=1:length(Data.trial)
    i
    v=Data.trial{1,i}(1,:);
    x=[x;v];
end

c1_veh=[];
for i=1:size(x,1)  
   if max(x(i,:))>48 || min(x(i,:))<-45
       c1_veh=[c1_veh;i];
   else
        plot(x(i,:),'b-')
        hold on
   end
        
end
ylim([-300 300])



GC_cluster1_veh_median_wa       =GC_cluster1_veh;
GC_Bp_cluster1_veh_median_wa    =GC_Bp_cluster1_veh;
GC_time_cluster1_veh_median_wa  =GC_time_cluster1_veh;
i= c1_veh';
GC_cluster1_veh_median_wa(i,:)       =[];
GC_Bp_cluster1_veh_median_wa(i,:)    =[];
GC_time_cluster1_veh_median_wa(i,:)  =[];


%save GC_cluster1_veh_median_wa.mat


%% RGS

%% per cluster
R = (cellfun(@(equis1) max(abs(hilbert(equis1(2,3001-50:3001+50)))),GC_Bp_cluster1_rgs));

[~,r_nl]=sort(abs(R-median(R)),'ascend');
R=R(r_nl);
GC_cluster1_rgs=GC_cluster1_rgs(r_nl);
GC_cluster1_rgs = GC_cluster1_rgs(1:2036);

GC_Bp_cluster1_rgs=GC_Bp_cluster1_rgs(r_nl);
GC_Bp_cluster1_rgs= GC_Bp_cluster1_rgs(1:2036);

GC_time_cluster1_rgs=GC_time_cluster1_rgs(r_nl,:);
GC_time_cluster1_rgs = GC_time_cluster1_rgs(1:2036,:);

fn=1000;
leng=length(GC_cluster1_rgs);
ro=3000;
tm = create_timecell(ro,leng);

label=[{'PFC'}; {'HPC'}];
Data.label=label;
Data.time=tm;
Data.trial=GC_cluster1_rgs.';

% NOTCH
cfg.bsfilter = 'yes';
cfg.bsfreq = [49 51];
Data = ft_preprocessing(cfg,Data);

% Run bandpass filter for plotting
cfg.bpfilter = 'yes';
cfg.bpfreq = [100 300];
Data = ft_preprocessing(cfg,Data);


% Detect artifact (overlay all ripples)
x=[];
for i=1:length(Data.trial)
    i
    v=Data.trial{1,i}(1,:);
    x=[x;v];
end

c1_rgs=[];
for i=1:size(x,1)  
   if max(x(i,:))>42 || min(x(i,:))<-48
       c1_rgs=[c1_rgs;i];
   else
        plot(x(i,:),'b-')
        hold on
   end      
end
ylim([-300 300])

GC_cluster4_rgs_median_wa       =GC_cluster1_rgs;
GC_Bp_cluster4_rgs_median_wa    =GC_Bp_cluster1_rgs;
GC_time_cluster4_rgs_median_wa  =GC_time_cluster1_rgs;
i=c1_rgs';
GC_cluster4_rgs_median_wa(i,:)       =[];
GC_Bp_cluster4_rgs_median_wa(i,:)    =[];
GC_time_cluster4_rgs_median_wa(i,:)  =[];



save('GC_ripple_4clusters_median_wa.mat','-v7.3')
