%% Checking number of pyr neurons match with excel sheet. 
% 101 RGS, 57 VEH. 

% RGS
cd('/home/adrian/Documents/GitHub/RGS14/spiking_activity')

load('SD_Data_RGS_Pyr_clean.mat')
NPD=[];
for i=1:length(SD_Data_RGS_Pyr) % Rats
    
    for ii=1: length(SD_Data_RGS_Pyr(i).RatxData)% Study days
        
         x=SD_Data_RGS_Pyr(i).RatxData(ii).SD_Spike_Info;
         xx=struct2table(x);
         xxx=table2cell(xx(:,1));
         
         NPD=[ NPD sum(~cellfun('isempty',xxx))]; % Neurons per day.
         
    end
    
end
sum(NPD)
%% VEH
load('SD_Data_Veh_Pyr_clean.mat')

NPD=[];
cont=0;
for i=1:length(SD_Data_Veh_Pyr) % Rats
    
    for ii=1: length(SD_Data_Veh_Pyr(i).Ratx)% Study days
        
[g2g.('Presleep'),Duration_off.('Presleep'),Count_off.('Presleep'),Rate_off.('Presleep'),Duration_on.('Presleep'),Count_on.('Presleep'),Rate_on.('Presleep'), NREMdur.('Presleep'),NPD]=get_off_data_trial(SD_Data_Veh_Pyr,i,ii,NPD,'Presleep');
[g2g.('Post_Trial_1'),Duration_off.('Post_Trial_1'),Count_off.('Post_Trial_1'),Rate_off.('Post_Trial_1'),Duration_on.('Post_Trial_1'),Count_on.('Post_Trial_1'),Rate_on.('Post_Trial_1'),NREMdur.('Post_Trial_1'),~]=get_off_data_trial(SD_Data_Veh_Pyr,i,ii,NPD,'Post_Trial_1');
[g2g.('Post_Trial_2'),Duration_off.('Post_Trial_2'),Count_off.('Post_Trial_2'),Rate_off.('Post_Trial_2'),Duration_on.('Post_Trial_2'),Count_on.('Post_Trial_2'),Rate_on.('Post_Trial_2'),NREMdur.('Post_Trial_2'),~]=get_off_data_trial(SD_Data_Veh_Pyr,i,ii,NPD,'Post_Trial_2');
[g2g.('Post_Trial_3'),Duration_off.('Post_Trial_3'),Count_off.('Post_Trial_3'),Rate_off.('Post_Trial_3'),Duration_on.('Post_Trial_3'),Count_on.('Post_Trial_3'),Rate_on.('Post_Trial_3'),NREMdur.('Post_Trial_3'),~]=get_off_data_trial(SD_Data_Veh_Pyr,i,ii,NPD,'Post_Trial_3');
[g2g.('Post_Trial_4'),Duration_off.('Post_Trial_4'),Count_off.('Post_Trial_4'),Rate_off.('Post_Trial_4'),Duration_on.('Post_Trial_4'),Count_on.('Post_Trial_4'),Rate_on.('Post_Trial_4'),NREMdur.('Post_Trial_4'),~]=get_off_data_trial(SD_Data_Veh_Pyr,i,ii,NPD,'Post_Trial_4');


[g2g.('Post_Trial_5_1'),Duration_off.('Post_Trial_5_1'),Count_off.('Post_Trial_5_1'),Rate_off.('Post_Trial_5_1'),Duration_on.('Post_Trial_5_1'),Count_on.('Post_Trial_5_1'),Rate_on.('Post_Trial_5_1'),NREMdur.('Post_Trial_5_1'),~]=get_off_data_trial_pt5(SD_Data_Veh_Pyr,i,ii,NPD,'PT5_Part_1');
[g2g.('Post_Trial_5_2'),Duration_off.('Post_Trial_5_2'),Count_off.('Post_Trial_5_2'),Rate_off.('Post_Trial_5_2'),Duration_on.('Post_Trial_5_2'),Count_on.('Post_Trial_5_2'),Rate_on.('Post_Trial_5_2'),NREMdur.('Post_Trial_5_2'),~]=get_off_data_trial_pt5(SD_Data_Veh_Pyr,i,ii,NPD,'PT5_Part_2');
[g2g.('Post_Trial_5_3'),Duration_off.('Post_Trial_5_3'),Count_off.('Post_Trial_5_3'),Rate_off.('Post_Trial_5_3'),Duration_on.('Post_Trial_5_3'),Count_on.('Post_Trial_5_3'),Rate_on.('Post_Trial_5_3'),NREMdur.('Post_Trial_5_3'),~]=get_off_data_trial_pt5(SD_Data_Veh_Pyr,i,ii,NPD,'PT5_Part_3');
[g2g.('Post_Trial_5_4'),Duration_off.('Post_Trial_5_4'),Count_off.('Post_Trial_5_4'),Rate_off.('Post_Trial_5_4'),Duration_on.('Post_Trial_5_4'),Count_on.('Post_Trial_5_4'),Rate_on.('Post_Trial_5_4'),NREMdur.('Post_Trial_5_4'),~]=get_off_data_trial_pt5(SD_Data_Veh_Pyr,i,ii,NPD,'PT5_Part_4');


cont=cont+1;
% 
% if cont==2
%     xo
% end
[Count_on_sum(cont,:)]=data_reduction(Count_on,'sum');
[Count_off_sum(cont,:)]=data_reduction(Count_off,'sum');

[NREMdur_sum(cont,:)]=data_reduction(NREMdur,'sum');

incidence_on(cont,:)=Count_on_sum(cont,:)./NREMdur_sum(cont,:);
incidence_off(cont,:)=Count_off_sum(cont,:)./NREMdur_sum(cont,:);

[Duration_on_mean(cont,:)]=data_reduction(Duration_on,'mean');
[Duration_off_mean(cont,:)]=data_reduction(Duration_off,'mean');

%xo
    end
    
end
sum(NPD)


[inde,ind]=sort(NPD,'descend')
%% Sort by neuron number
NPD=NPD(ind)
incidence_on=incidence_on(ind,:);
incidence_off=incidence_off(ind,:);
Duration_on_mean=Duration_on_mean(ind,:);
Duration_off_mean=Duration_off_mean(ind,:);
%%
Duration_on_mean=Duration_on_mean.*1000; %ms
Duration_off_mean=Duration_off_mean.*1000;
%%
% RGS
cd('/home/adrian/Documents/GitHub/RGS14/spiking_activity')

load('SD_Data_RGS_Pyr_clean.mat')

NPD=[];
cont=0;
for i=1:length(SD_Data_RGS_Pyr) % Rats
    
    for ii=1: length(SD_Data_RGS_Pyr(i).RatxData)% Study days
        
[g2g.('Presleep'),Duration_off.('Presleep'),Count_off.('Presleep'),Rate_off.('Presleep'),Duration_on.('Presleep'),Count_on.('Presleep'),Rate_on.('Presleep'), NREMdur.('Presleep'),NPD]=get_off_data_trial(SD_Data_RGS_Pyr,i,ii,NPD,'Presleep');
[g2g.('Post_Trial_1'),Duration_off.('Post_Trial_1'),Count_off.('Post_Trial_1'),Rate_off.('Post_Trial_1'),Duration_on.('Post_Trial_1'),Count_on.('Post_Trial_1'),Rate_on.('Post_Trial_1'),NREMdur.('Post_Trial_1'),~]=get_off_data_trial(SD_Data_RGS_Pyr,i,ii,NPD,'Post_Trial_1');
[g2g.('Post_Trial_2'),Duration_off.('Post_Trial_2'),Count_off.('Post_Trial_2'),Rate_off.('Post_Trial_2'),Duration_on.('Post_Trial_2'),Count_on.('Post_Trial_2'),Rate_on.('Post_Trial_2'),NREMdur.('Post_Trial_2'),~]=get_off_data_trial(SD_Data_RGS_Pyr,i,ii,NPD,'Post_Trial_2');
[g2g.('Post_Trial_3'),Duration_off.('Post_Trial_3'),Count_off.('Post_Trial_3'),Rate_off.('Post_Trial_3'),Duration_on.('Post_Trial_3'),Count_on.('Post_Trial_3'),Rate_on.('Post_Trial_3'),NREMdur.('Post_Trial_3'),~]=get_off_data_trial(SD_Data_RGS_Pyr,i,ii,NPD,'Post_Trial_3');
[g2g.('Post_Trial_4'),Duration_off.('Post_Trial_4'),Count_off.('Post_Trial_4'),Rate_off.('Post_Trial_4'),Duration_on.('Post_Trial_4'),Count_on.('Post_Trial_4'),Rate_on.('Post_Trial_4'),NREMdur.('Post_Trial_4'),~]=get_off_data_trial(SD_Data_RGS_Pyr,i,ii,NPD,'Post_Trial_4');


[g2g.('Post_Trial_5_1'),Duration_off.('Post_Trial_5_1'),Count_off.('Post_Trial_5_1'),Rate_off.('Post_Trial_5_1'),Duration_on.('Post_Trial_5_1'),Count_on.('Post_Trial_5_1'),Rate_on.('Post_Trial_5_1'),NREMdur.('Post_Trial_5_1'),~]=get_off_data_trial_pt5(SD_Data_RGS_Pyr,i,ii,NPD,'PT5_Part_1');
[g2g.('Post_Trial_5_2'),Duration_off.('Post_Trial_5_2'),Count_off.('Post_Trial_5_2'),Rate_off.('Post_Trial_5_2'),Duration_on.('Post_Trial_5_2'),Count_on.('Post_Trial_5_2'),Rate_on.('Post_Trial_5_2'),NREMdur.('Post_Trial_5_2'),~]=get_off_data_trial_pt5(SD_Data_RGS_Pyr,i,ii,NPD,'PT5_Part_2');
[g2g.('Post_Trial_5_3'),Duration_off.('Post_Trial_5_3'),Count_off.('Post_Trial_5_3'),Rate_off.('Post_Trial_5_3'),Duration_on.('Post_Trial_5_3'),Count_on.('Post_Trial_5_3'),Rate_on.('Post_Trial_5_3'),NREMdur.('Post_Trial_5_3'),~]=get_off_data_trial_pt5(SD_Data_RGS_Pyr,i,ii,NPD,'PT5_Part_3');
[g2g.('Post_Trial_5_4'),Duration_off.('Post_Trial_5_4'),Count_off.('Post_Trial_5_4'),Rate_off.('Post_Trial_5_4'),Duration_on.('Post_Trial_5_4'),Count_on.('Post_Trial_5_4'),Rate_on.('Post_Trial_5_4'),NREMdur.('Post_Trial_5_4'),~]=get_off_data_trial_pt5(SD_Data_RGS_Pyr,i,ii,NPD,'PT5_Part_4');


cont=cont+1;
% 
% if cont==2
%     xo
% end
[Count_on_sum(cont,:)]=data_reduction(Count_on,'sum');
[Count_off_sum(cont,:)]=data_reduction(Count_off,'sum');

[NREMdur_sum(cont,:)]=data_reduction(NREMdur,'sum');

incidence_on(cont,:)=Count_on_sum(cont,:)./NREMdur_sum(cont,:);
incidence_off(cont,:)=Count_off_sum(cont,:)./NREMdur_sum(cont,:);

[Duration_on_mean(cont,:)]=data_reduction(Duration_on,'mean');
[Duration_off_mean(cont,:)]=data_reduction(Duration_off,'mean');

%xo
    end
    
end
sum(NPD)


 [inde,ind]=sort(NPD,'descend')
%% Sort by neuron number
NPD=NPD(ind)
incidence_on=incidence_on(ind,:);
incidence_off=incidence_off(ind,:);
Duration_on_mean=Duration_on_mean(ind,:);
Duration_off_mean=Duration_off_mean(ind,:);
%%
Duration_on_mean=Duration_on_mean.*1000; %ms
Duration_off_mean=Duration_off_mean.*1000;

