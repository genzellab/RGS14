%% MAIN
clc
clear
format compact
format longG
load('waveforms_ripple_clusters_all.mat')

%save('waveforms_cluster2_rgs.mat','waveforms_cluster2_bp_rgs','waveforms_cluster2_raw_rgs','-v7.3');
%save('waveforms_cluster2_rgs.mat','waveforms_cluster2_bp_rgs','waveforms_cluster2_raw_rgs','-v7.3');
%save('waveforms_cluster2_rgs.mat','waveforms_cluster2_bp_rgs','waveforms_cluster2_raw_rgs','-v7.3');
%save('waveforms_cluster2_rgs.mat','waveforms_cluster2_bp_rgs','waveforms_cluster2_raw_rgs','-v7.3');


%% bandpassed
data=waveforms_cluster2_bp_rgs;
num_rat=unique(cell2mat(data(:,6)));
num_SD=unique(cell2mat(data(:,7)));
% RUN UNTIL HERE AT FIRST

% ACCORDING TO EACH FILE, ADJUST RATS AND SDs MANUALLY
rat3=cell(1,10);
rat4=cell(1,10);
rat7=cell(1,10);
rat8=cell(1,10);
for i=1:size(data,1)
    % rat number
    i
    if cell2mat(data(i,6))==3
        rat3_=data(i,:);
        rat3=[rat3;rat3_];
    elseif cell2mat(data(i,6))==4
        rat4_=data(i,:);
        rat4=[rat4;rat4_];
    elseif cell2mat(data(i,6))==7
        rat7_=data(i,:);
        rat7=[rat7;rat7_];
    elseif cell2mat(data(i,6))==8
        rat8_=data(i,:);
        rat8=[rat8;rat8_];
    end
end



%%3 rat
x=(cell2mat(rat3(:,7))==1);
rat3_SD1_bp_c2_rgs=rat3(x==1,:);

x=(cell2mat(rat3(:,7))==2);
rat3_SD2_bp_c2_rgs=rat3(x==1,:);

x=(cell2mat(rat3(:,7))==3);
rat3_SD3_bp_c2_rgs=rat3(x==1,:);

x=(cell2mat(rat3(:,7))==4);
rat3_SD4_bp_c2_rgs=rat3(x==1,:);

x=(cell2mat(rat3(:,7))==5);
rat3_SD5_bp_c2_rgs=rat3(x==1,:);

x=(cell2mat(rat3(:,7))==6);
rat3_SD6_bp_c2_rgs=rat3(x==1,:);

x=(cell2mat(rat3(:,7))==14);
rat3_SD14_bp_c2_rgs=rat3(x==1,:);



%%4 rat
x=(cell2mat(rat4(:,7))==1);
rat4_SD1_bp_c2_rgs=rat4(x==1,:);

x=(cell2mat(rat4(:,7))==2);
rat4_SD2_bp_c2_rgs=rat4(x==1,:);

x=(cell2mat(rat4(:,7))==3);
rat4_SD3_bp_c2_rgs=rat4(x==1,:);

x=(cell2mat(rat4(:,7))==4);
rat4_SD4_bp_c2_rgs=rat4(x==1,:);

x=(cell2mat(rat4(:,7))==5);
rat4_SD5_bp_c2_rgs=rat4(x==1,:);

x=(cell2mat(rat4(:,7))==6);
rat4_SD6_bp_c2_rgs=rat4(x==1,:);

x=(cell2mat(rat4(:,7))==14);
rat4_SD14_bp_c2_rgs=rat4(x==1,:);




%%7 rat
x=(cell2mat(rat7(:,7))==1);
rat7_SD1_bp_c2_rgs=rat7(x==1,:);

x=(cell2mat(rat7(:,7))==2);
rat7_SD2_bp_c2_rgs=rat7(x==1,:);

x=(cell2mat(rat7(:,7))==3);
rat7_SD3_bp_c2_rgs=rat7(x==1,:);

x=(cell2mat(rat7(:,7))==4);
rat7_SD4_bp_c2_rgs=rat7(x==1,:);

x=(cell2mat(rat7(:,7))==5);
rat7_SD5_bp_c2_rgs=rat7(x==1,:);

x=(cell2mat(rat7(:,7))==6);
rat7_SD6_bp_c2_rgs=rat7(x==1,:);

x=(cell2mat(rat7(:,7))==14);
rat7_SD14_bp_c2_rgs=rat7(x==1,:);




%%8 rat
x=(cell2mat(rat8(:,7))==1);
rat8_SD1_bp_c2_rgs=rat8(x==1,:);

x=(cell2mat(rat8(:,7))==2);
rat8_SD2_bp_c2_rgs=rat8(x==1,:);

x=(cell2mat(rat8(:,7))==3);
rat8_SD3_bp_c2_rgs=rat8(x==1,:);

x=(cell2mat(rat8(:,7))==4);
rat8_SD4_bp_c2_rgs=rat8(x==1,:);

x=(cell2mat(rat8(:,7))==5);
rat8_SD5_bp_c2_rgs=rat8(x==1,:);

x=(cell2mat(rat8(:,7))==6);
rat8_SD6_bp_c2_rgs=rat8(x==1,:);

x=(cell2mat(rat8(:,7))==14);
rat8_SD14_bp_c2_rgs=rat8(x==1,:);


save('rat_SD_grouped_cluster2_rgs.mat');


%%find empty cells
pp=cellfun(@isempty,GC_raw_veh(:,7));
[i,j] = find(pp);

% replace with a number
start=51938;
end_=52151;
sd=2;
GC_raw_veh(start:end_,7)=num2cell(sd);
GC_bandpassed_veh(start:end_,7)=num2cell(sd);

%
save('waveforms_ripple_cluster2_rgs.mat','waveforms_cluster2_bp_rgs','waveforms_cluster2_raw_rgs','-v7.3');



