function [granger,granger2]=granger_automation(input1,input2,freqrange)

label = [{'PFC'}; {'HPC'}];

%%  input1
clear Data
fn=1000;
leng=length(input1);
ro=3000;
tm=create_timecell(ro,leng);

Data.label=label;
Data.time=tm;
Data.trial=input1.';
%
cfg.bsfilter = 'yes';
cfg.bsfreq = [49 51];
Data = ft_preprocessing(cfg,Data);

[granger,freq]=createauto_timefreq(Data,freqrange,[-1.1:0.01:1.1]);

%%  input2
clear Data
fn = 1000;
leng = length(input2);
ro = 3000;
tm = create_timecell(ro,leng);

Data.label = label;
Data.time = tm;
Data.trial = input2.';
%
cfg.bsfilter = 'yes';
cfg.bsfreq = [49 51];
Data = ft_preprocessing(cfg,Data);

[granger2,freq2]=createauto_timefreq(Data,freqrange,[-1.1:0.01:1.1]);

end