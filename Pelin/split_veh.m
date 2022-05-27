%% Rat 1 

Rat_i = find([waveforms_cluster1_raw_veh{:,6}]==1);
Rat1_raw_veh = waveforms_cluster1_raw_veh(Rat_i,:);
%save ('Rat1_compiled_c1','Rat1_raw_veh')

SD_i = find([Rat1_raw_veh{:,7}]==1);
Rat1_SD1_c1_raw_veh = Rat1_raw_veh(SD_i,:);
save('Rat1_SD1_CON_c1_raw_veh','Rat1_SD1_c1_raw_veh')

SD_i = find([Rat1_raw_veh{:,7}]==2);
Rat1_SD2_c1_raw_veh = Rat1_raw_veh(SD_i,:);
save('Rat1_SD2_OD_c1_raw_veh','Rat1_SD2_c1_raw_veh')

SD_i = find([Rat1_raw_veh{:,7}]==3);
Rat1_SD3_c1_raw_veh = Rat1_raw_veh(SD_i,:);
save('Rat1_SD3_OR_c1_raw_veh','Rat1_SD3_c1_raw_veh')

SD_i = find([Rat1_raw_veh{:,7}]==4);
Rat1_SD4_c1_raw_veh = Rat1_raw_veh(SD_i,:);
save('Rat1_SD4_HC_c1_raw_veh','Rat1_SD4_c1_raw_veh')

%% Rat 2

Rat_i = find([waveforms_cluster1_raw_veh{:,6}]==2);
Rat2_raw_veh = waveforms_cluster1_raw_veh(Rat_i,:);
%save ('Rat2_compiled_c1_raw_veh','Rat2_raw_veh')

SD_i = find([Rat2_raw_veh{:,7}]==1);
Rat2_SD1_c1_raw_veh = Rat2_raw_veh(SD_i,:);
save('Rat2_SD1_OD_c1_raw_veh','Rat2_SD1_c1_raw_veh')

SD_i = find([Rat2_raw_veh{:,7}]==2);
Rat2_SD2_c1_raw_veh = Rat2_raw_veh(SD_i,:);
save('Rat2_SD2_OR_c1_raw_veh','Rat2_SD2_c1_raw_veh')

SD_i = find([Rat2_raw_veh{:,7}]==3);
Rat2_SD3_c1_raw_veh = Rat2_raw_veh(SD_i,:);
save('Rat2_SD3_CON_c1_raw_veh','Rat2_SD3_c1_raw_veh')

SD_i = find([Rat2_raw_veh{:,7}]==4);
Rat2_SD4_c1_raw_veh = Rat2_raw_veh(SD_i,:);
save('Rat2_SD4_HC_c1_raw_veh','Rat2_SD4_c1_raw_veh')

%% Rat 6

Rat_i = find([waveforms_cluster1_raw_veh{:,6}]==6);
Rat6_raw_veh = waveforms_cluster1_raw_veh(Rat_i,:);

SD_i = find(contains(Rat6_raw_veh(:,8),'HC'));
Rat6_raw_veh(:,7) =  num2cell(1);
Rat6_SD1_c1_raw_veh = Rat6_raw_veh(SD_i,:);
save('Rat6_SD1_HC_c1_raw_veh','Rat6_SD1_c1_raw_veh')

SD_i = find(contains(Rat6_raw_veh(:,8),'OR'));
Rat6_raw_veh(:,7) =  num2cell(2);
Rat6_SD2_c1_raw_veh = Rat6_raw_veh(SD_i,:);
save('Rat6_SD2_OR_c1_raw_veh','Rat6_SD2_c1_raw_veh')

SD_i = find(contains(Rat6_raw_veh(:,8),'CON'));
Rat6_raw_veh(:,7) =  num2cell(3);
Rat6_SD3_c1_raw_veh = Rat6_raw_veh(SD_i,:);
save('Rat6_SD3_CON_c1_raw_veh','Rat6_SD3_c1_raw_veh')

SD_i = find(contains(Rat6_raw_veh(:,8),'OD'));
Rat6_raw_veh(:,7) =  num2cell(4);
Rat6_SD4_c1_raw_veh = Rat6_raw_veh(SD_i,:);
save('Rat6_SD4_OD_c1_raw_veh','Rat6_SD4_c1_raw_veh')

%save('rat6_compiled_c1_raw_veh','Rat6_raw_veh')


%% Rat 9

Rat_i = find([waveforms_cluster1_raw_veh{:,6}]==9);
Rat9_raw_veh = waveforms_cluster1_raw_veh(Rat_i,:);
%save('Rat9_compiled_c1_raw_veh','Rat9_raw_veh')

SD_i = find([Rat9_raw_veh{:,7}]==1);
Rat9_SD1_c1_raw_veh = Rat9_raw_veh(SD_i,:);
save('Rat9_SD1_HC_c1_raw_veh','Rat9_SD1_c1_raw_veh')

SD_i = find([Rat9_raw_veh{:,7}]==2);
Rat9_SD2_c1_raw_veh = Rat9_raw_veh(SD_i,:);
save('Rat9_SD2_CON_c1_raw_veh','Rat9_SD2_c1_raw_veh')

SD_i = find([Rat9_raw_veh{:,7}]==3);
Rat9_SD3_c1_raw_veh = Rat9_raw_veh(SD_i,:);
save('Rat9_SD3_OR_c1_raw_veh','Rat9_SD3_c1_raw_veh')

SD_i = find([Rat9_raw_veh{:,7}]==6);
Rat9_SD6_c1_raw_veh = Rat9_raw_veh(SD_i,:);
save('Rat9_SD6_OD_c1_raw_veh','Rat9_SD6_c1_raw_veh')
