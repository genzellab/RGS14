function [g2g,Duration_off,Count_off,Rate_off,Duration_on,Count_on,Rate_on,NonREMdur_tot,NPD]=get_off_data_trial(SD_Data_Veh_Pyr,i,ii,NPD,Trial)

         x=SD_Data_Veh_Pyr(i).RatxData(ii).SD_Spike_Info;
         xx=struct2table(x);
         xxx=table2cell(xx(:,1));
         
         NPD=[ NPD sum(~cellfun('isempty',xxx))]; % Neurons per day
         
 Count_off=[];
 Duration_off=[];
 Rate_off=[];
 
 Shuffled=[];
 
         
 Count_on=[];
 Duration_on=[];
 Rate_on=[];
 NonREMdur_tot=[];
 
 jk=1;
 if isempty(x(jk).(Trial))
     jk=2;
 end
 uplim=length(x(jk).(Trial).NREM_Spikes);
 %uplim
         parfor iii=1:uplim % NREM bouts. 
             X=[]; % compressed data (mua) per bout
             Y=[];
             for iiii=1:length(xxx) % rows/neurons. (Some rows are empty)

                if isempty(x(iiii).(Trial))
                    continue
                end
                X=[X; x(iiii).(Trial).NREM_Spikes{iii,1}];
                y1=x(iiii).(Trial).NREM_Bout_Info{iii,1};
                y2=x(iiii).(Trial).NREM_Bout_Info{iii,2};
                Y=[Y; str2num(y1);str2num(y2)];

             end
             X=sort(unique(X));
             Y=unique(Y);
             if length(Y)~=2
                error('Multiple start and end timesstamp for bout. ')
             end
%              stem(X, (X*0)+50);hold on
%               stem(Y,(Y*0)+70)
              BI=ismember(Y(1):Y(2),X);
%%
            %on-off  
[count_off,duration_off,rate_off,count_on,duration_on,rate_on,NonREMdur]=find_off_period(BI,30000,0.050);            

if strcmp(Trial,'Presleep1')
    Mean_dur_shuf=[];
    for shuf=1:500
    [~,duration_off_shuf]=find_off_period(BI(randperm(length(BI))),30000,0.050);          
    Mean_dur_shuf=[Mean_dur_shuf mean(duration_off_shuf)];

    end
    Shuffled=[Shuffled;Mean_dur_shuf ];
end

Count_off=[Count_off; count_off ];
Duration_off=[Duration_off duration_off];
Rate_off=[Rate_off rate_off];

Count_on=[Count_on; count_on ];
Duration_on=[Duration_on duration_on];
Rate_on=[Rate_on rate_on];

NonREMdur_tot=[NonREMdur_tot NonREMdur];
%%
             
             
% histogram(duration_on*1000,[50:10:1000])
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
% xlim([50 1000])
% 
% hold on
% histogram(duration_on_shuf*1000,[50:10:1000])
% set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')
%%
% mean(duration_on*1000)
% up_bound = prctile(Mean_dur_shuf*1000,97.5)
%%
             %%
         end
if ~isempty(Shuffled)
    if mean(Duration_off*1000)>  prctile(Shuffled(:)*1000,97.5) 
    g2g=1; %Good to go.
    else
    g2g=0;
    error('Not above chance')

    end
else
    g2g=[];
end
  %Duration_off=mean(Duration_off);
  %Count_off=sum(Count_off);
end