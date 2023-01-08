function [freq1,freq2,subs_freq,zlim]=spectrogram_automation(input1,input2,channel,freqrange)

bottom=min([size(input1,1) size(input2,1)]);
input1=input1(randperm(length(input1)));
input1=input1(1:bottom);
input2=input2(randperm(length(input2)));
input2=input2(1:bottom);

%% input1
clear Data
fn=1000;
leng=length(input1);
ro=3000;
tm = create_timecell(ro,leng);

label=[{'PFC'}; {'HPC'}];
Data.label=label;
Data.time=tm;
Data.trial=input1.';

% Notch filter
cfg.bsfilter = 'yes';
cfg.bsfreq = [49 51];
Data = ft_preprocessing(cfg,Data);

[freq1]=time_frequency(Data,freqrange,[-1.1:0.01:1.1]);  % use 10 ms for the analysis

%% input2
clear Data
fn=1000;
leng=length(input2);
ro=3000;
tm = create_timecell(ro,leng);

label=[{'PFC'}; {'HPC'}];
Data.label=label;
Data.time=tm;
Data.trial=input2.';

% Notch filter
cfg.bsfilter = 'yes';
cfg.bsfreq = [49 51];
Data = ft_preprocessing(cfg,Data);

[freq2]=time_frequency(Data,freqrange,[-1.1:0.01:1.1]);  % use 10 ms for the analysis

%% Contrast
cfg_m = [];
cfg_m.operation = 'subtract';
cfg_m.parameter = 'powspctrm';
subs_freq = ft_math(cfg_m,  freq1, freq2);

%%
cfg              = [];
cfg.channel      = channel;
[ zmin1, zmax1] = ft_getminmax(cfg, freq1);
[zmin2, zmax2] = ft_getminmax(cfg, freq2);
zlim=[min([zmin1 zmin2]) max([zmax1 zmax2])];

end