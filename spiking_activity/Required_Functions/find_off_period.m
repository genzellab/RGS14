function [count_off,duration_off,rate_off,count_on,duration_on,rate_on,NonREMdur]=find_off_period(BI,fn,duration_min)
%BI: Binary vector representing continuous timestamps with 1 meaning a
%spike was detected for this timestamp. 
% fn: sampling rate
            NonREMdur=length(BI)/fn; %sec

            [spikes]=myFind(BI==1);
            spikes=spikes.';
            isi=diff(spikes)/fn; % Inter spike interval (sec)
            off_ind=myFind(isi>duration_min);
            off_ind=off_ind.';
            
            off_start=spikes(off_ind);
            off_end=spikes(off_ind+1);

            count_off=length(off_start);
            duration_off=isi(off_ind);
            rate_off=count_off/NonREMdur;

        %ON period
        on_s=[off_start(2:end);off_end(1:end-1)];
        on_start=on_s(2,[off_start(2:end)-off_end(1:end-1)]~=0);
        on_end=on_s(1,[off_start(2:end)-off_end(1:end-1)]~=0);
        on_start_new=on_start(myFind([(on_end-on_start)/30]>=50 & [(on_end-on_start)/30]<=4000 ));
        on_end_new=on_end(myFind([(on_end-on_start)/30]>=50 & [(on_end-on_start)/30]<=4000));
        
        left_spikes=spikes(~ismember(spikes,off_start) & ~ismember(spikes,off_end) & ~ismember(spikes,on_end_new) & ~ismember(spikes,on_start_new));
        
        [n_on_start_new, n_on_end_new]=on_minimum_spikenum(on_start_new, on_end_new,left_spikes);

        count_on=length(n_on_start_new);
        duration_on=(n_on_end_new-n_on_start_new)/fn;
            if isempty(duration_on)
                duration_on=[];
            end
        rate_on=count_on/NonREMdur;

end