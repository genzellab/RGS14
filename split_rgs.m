%% Rat 3 
Rat_i = find([waveforms_cluster1_raw_rgs{:,6}]==3);
Rat3_raw_rgs = waveforms_cluster1_raw_rgs(Rat_i,:);
%save('Rat3_compiled_c1_raw_rgs.mat','Rat3_raw_rgs')

SD_i = find([Rat3_raw_rgs{:,7}]==1);
Rat3_SD1_c1_raw_rgs = Rat3_raw_rgs(SD_i,:);
save('Rat3_SD1_OD_c1_raw_rgs.mat','Rat3_SD1_c1_raw_rgs')

SD_i = find([Rat3_raw_rgs{:,7}]==3);
Rat3_SD3_c1_raw_rgs = Rat3_raw_rgs(SD_i,:);
save('Rat3_SD3_OR_c1_raw_rgs.mat','Rat3_SD3_c1_raw_rgs')

SD_i = find([Rat3_raw_rgs{:,7}]==5);
Rat3_SD5_c1_raw_rgs = Rat3_raw_rgs(SD_i,:);
save('Rat3_SD5_CON_c1_raw_rgs.mat','Rat3_SD5_c1_raw_rgs')

SD_i = find([Rat3_raw_rgs{:,7}]==14);
Rat3_SD14_c1_raw_rgs = Rat3_raw_rgs(SD_i,:);
save('Rat3_SD14_HC_c1_raw_rgs.mat','Rat3_SD14_c1_raw_rgs')

%% Rat 4
Rat_i = find([waveforms_cluster1_raw_rgs{:,6}]==4);
Rat4_raw_rgs = waveforms_cluster1_raw_rgs(Rat_i,:);
%save('Rat4_compiled_c1_raw_rgs.mat','Rat4_raw_rgs')

SD_i = find([Rat4_raw_rgs{:,7}]==3);
Rat4_SD3_c1_raw_rgs = Rat4_raw_rgs(SD_i,:);
save('Rat4_SD3_CON_c1_raw_rgs.mat','Rat4_SD3_c1_raw_rgs')

SD_i = find([Rat4_raw_rgs{:,7}]==4);
Rat4_SD4_c1_raw_rgs = Rat4_raw_rgs(SD_i,:);
save('Rat4_SD4_HC_c1_raw_rgs.mat','Rat4_SD4_c1_raw_rgs')

SD_i = find([Rat4_raw_rgs{:,7}]==5);
Rat4_SD5_c1_raw_rgs = Rat4_raw_rgs(SD_i,:);
save('Rat4_SD5_OD_c1_raw_rgs.mat','Rat4_SD5_c1_raw_rgs')

SD_i = find([Rat4_raw_rgs{:,7}]==6);
Rat4_SD6_c1_raw_rgs = Rat4_raw_rgs(SD_i,:);
save('Rat4_SD6_OR_c1_raw_rgs.mat','Rat4_SD6_c1_raw_rgs')

%% Rat 7
Rat_i = find([waveforms_cluster1_raw_rgs{:,6}]==7);
Rat7_raw_rgs = waveforms_cluster1_raw_rgs(Rat_i,:);
%save('Rat7_compiled_c1_raw_rgs.mat','Rat7_raw_rgs')

SD_i = find([Rat7_raw_rgs{:,7}]==1);
Rat7_SD1_c1_raw_rgs = Rat7_raw_rgs(SD_i,:);
save('Rat7_SD1_HC_c1_raw_rgs.mat','Rat7_SD1_c1_raw_rgs')

SD_i = find([Rat7_raw_rgs{:,7}]==2);
Rat7_SD2_c1_raw_rgs = Rat7_raw_rgs(SD_i,:);
save('Rat7_SD2_OR_c1_raw_rgs.mat','Rat7_SD2_c1_raw_rgs')

SD_i = find([Rat7_raw_rgs{:,7}]==3);
Rat7_SD3_c1_raw_rgs = Rat7_raw_rgs(SD_i,:);
save('Rat7_SD3_CON_c1_raw_rgs.mat','Rat7_SD3_c1_raw_rgs')

SD_i = find([Rat7_raw_rgs{:,7}]==5);
Rat7_SD5_c1_raw_rgs = Rat7_raw_rgs(SD_i,:);
save('Rat7_SD5_OD_c1_raw_rgs.mat','Rat7_SD5_c1_raw_rgs')

%% Rat 8 

Rat_i = find([waveforms_cluster1_raw_rgs{:,6}]==8);
Rat8_raw_rgs = waveforms_cluster1_raw_rgs(Rat_i,:);
%save('Rat8_compiled_c1_raw_rgs.mat','Rat8_raw_rgs')

SD_i = find([Rat8_raw_rgs{:,7}]==1);
Rat8_SD1_c1_raw_rgs = Rat8_raw_rgs(SD_i,:);
save('Rat8_SD1_HC_c1_raw_rgs.mat','Rat8_SD1_c1_raw_rgs')

SD_i = find([Rat8_raw_rgs{:,7}]==2);
Rat8_SD2_c1_raw_rgs = Rat8_raw_rgs(SD_i,:);
save('Rat8_SD2_CON_c1_raw_rgs.mat','Rat8_SD2_c1_raw_rgs')

SD_i = find([Rat8_raw_rgs{:,7}]==3);
Rat8_SD3_c1_raw_rgs = Rat8_raw_rgs(SD_i,:);
save('Rat8_SD3_OR_c1_raw_rgs.mat','Rat8_SD3_c1_raw_rgs')

SD_i = find([Rat8_raw_rgs{:,7}]==6);
Rat8_SD6_c1_raw_rgs = Rat8_raw_rgs(SD_i,:);
save('Rat8_SD6_OD_c1_raw_rgs.mat','Rat8_SD6_c1_raw_rgs')


