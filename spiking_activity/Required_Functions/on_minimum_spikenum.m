function [n_on_start_new, n_on_end_new]=on_minimum_spikenum(on_start_new, on_end_new,left_spikes)

n_on_start_new=on_start_new;
n_on_end_new=on_end_new;

%%
for u=1:length(on_start_new)
       interval=[ on_start_new(u) on_end_new(u) ];
      ai = myFind(left_spikes>=interval(1)&left_spikes<=interval(2));
%         length(ai)

        if length(ai)<8 % Including start and end would be 10 spikes. 
            n_on_start_new(u)=NaN;
            n_on_end_new(u)=NaN;
            
        end
end

n_on_start_new=n_on_start_new(~isnan(n_on_start_new));
n_on_end_new=n_on_end_new(~isnan(n_on_end_new));


%%
end