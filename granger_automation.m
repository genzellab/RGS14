function [granger,granger2,add]=granger_automation(input1,input2,input3,freqrange)

    label = [{'PFC'}; {'HPC'}];
    
    %%  input1
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
    
    %
    freqrange=[0:0.5:20];
    
    [granger,freq]=createauto_timefreq(Data,freqrange,[-1.1:0.01:1.1]);
    
    
    %%  input2
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
    
    %
    freqrange=[0:0.5:20];
    
    [granger2,freq2]=createauto_timefreq(Data,freqrange,[-1.1:0.01:1.1]);
    
    %% Contrast
    fn = 1000;
    leng = length(input3);
    ro = 3000;
    tm = create_timecell(ro,leng);
    
    Data.label = label;
    Data.time = tm;
    Data.trial = input3.';
    %
    cfg.bsfilter = 'yes';
    cfg.bsfreq = [49 51];
    Data = ft_preprocessing(cfg,Data); 
    
    %
    freqrange=[0:0.5:20];
    
    [add,freq3]=createauto_timefreq(Data,freqrange,[-1.1:0.01:1.1]);

end