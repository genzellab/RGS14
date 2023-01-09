function [granger,freq]=createauto_timefreq(data1,freqrange,toy)

equis=0.5;
       
cfg = [];
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.pad=8;
cfg.foi = freqrange;
cfg.t_ftimwin = 0.5.*ones(size(4./cfg.foi')); 
cfg.tapsmofrq = equis*cfg.foi;
cfg.toi=toy;
cfg.output         = 'fourier';
cfg.keeptrials = 'yes';
%
freq = ft_freqanalysis(cfg, data1);
%
cfg = [];
cfg.method    = 'granger';
granger    = ft_connectivityanalysis(cfg, freq);

end
