clear variables
addpath '/Volumes/Samsung_T5/Milan_DA/RGS14_Ephys_da/Data_RGS14_Downsampled_First_Session'
cd('/Volumes/Samsung_T5/Milan_DA/RGS14_Ephys_da/Data_RGS14_Downsampled_First_Session')
%% Get the study day folders
rat_folder=getfolder;
prompt = {'Enter the rat index'};
dlgtitle = 'Rat Index';
k = str2double(inputdlg(prompt,dlgtitle));
cd(rat_folder{k}) 
SD_folders  = getfolder;

%% Loading the timestamps and initiating variables
for j = 1:length(SD_folders)
    
    load(strcat('ripple_total_data_',SD_folders{j},'.mat'))
    load(strcat('spindles_total_data_',SD_folders{j},'.mat'))
    load(strcat('delta_timestamps_',SD_folders{j},'.mat'))
    delta_total_data = delta_timestamps_SD;
    
        current_sd = SD_folders{j};
        current_sd = strsplit(current_sd, '_');
        idx = find(contains(current_sd, 'sd', 'IgnoreCase',true));
        SDn = current_sd{idx}(1,3:end);
        cd ('/Volumes/Samsung_T5/Milan_DA/SWR_Cluster_SDwise_pelin')
        ClSDw =  getfolder;
        cd(ClSDw{k})
        path = cd;
        dinfo = dir(path);
        dinfo = {dinfo.name};
        dinfo = dinfo(4:end);
        idx2 = find(contains(dinfo,'trialwise_', 'IgnoreCase',true));
        dinfo = dinfo(idx2);
        idx3 = find(contains(dinfo, ['SD',SDn,'_'], 'IgnoreCase',true));
        dinfo_sd = dinfo(idx3);
       clus1 = load(dinfo_sd{contains(dinfo_sd,'c1')});
       name = fieldnames(clus1);
       clus1 = clus1.(name{2});
       clus2 = load(dinfo_sd{contains(dinfo_sd,'c2')});
       name = fieldnames(clus2);
       clus2 = clus2.(name{2});
       clus3 = load(dinfo_sd{contains(dinfo_sd,'c3')});
        name = fieldnames(clus3);
       clus3 = clus3.(name{2});
       clus4 = load(dinfo_sd{contains(dinfo_sd,'c4')});
        name = fieldnames(clus4);
       clus4 = clus4.(name{2});  
       
cd('/Volumes/Samsung_T5/Milan_DA/RGS14_Ephys_da/Data_RGS14_Downsampled_First_Session') 

cd(rat_folder{k})  
CO_ripples_c1 = {};
CO_spindles_c1 = {};
CO_count_hpc_c1 = [];
CO_count_pfc_c1 = [];

CO_ripples_c2 = {};
CO_spindles_c2 = {};
CO_count_hpc_c2 = [];
CO_count_pfc_c2 = [];

CO_ripples_c3 = {};
CO_spindles_c3 = {};
CO_count_hpc_c3 = [];
CO_count_pfc_c3 = [];

CO_ripples_c4 = {};
CO_spindles_c4 = {};
CO_count_hpc_c4 = [];
CO_count_pfc_c4 = [];

%% Extraction of indices and counts wrt Co-occurring events
    for i = 1:length(spindles_total_data)
        if ~isnan(spindles_total_data{1,i})
        spindles = spindles_total_data{i};
        deltas = delta_total_data{i};    
        ripples = ripple_total_data{i};
        ripples_c1 = cell2mat(clus1{i}(:,2:4));
        ripples_c2 = cell2mat(clus2{i}(:,2:4));
        ripples_c3 = cell2mat(clus3{i}(:,2:4));
        ripples_c4 = cell2mat(clus4{i}(:,2:4));
  
 %% Extracting Concatenated NREM timestamps per cluster 
       
if ~isnan(ripples)
    
        if ~isempty(ripples_c1)
        idx4 = ismember(ripples(:,1:3),ripples_c1,'rows'); 
        con_r_c1 = ripples(idx4, 5:7);                   % Cluster 1
        else
         con_r_c1 = [];   
        end
        
        if ~isempty(ripples_c2)
        idx5 = ismember(ripples(:,1:3),ripples_c2,'rows');
        con_r_c2 = ripples(idx5, 5:7);                   % Cluster 2
        else
         con_r_c2 = [];   
        end
        
        if ~isempty(ripples_c3)
        idx6 = ismember(ripples(:,1:3),ripples_c3,'rows');
        con_r_c3 = ripples(idx6, 5:7);                   % Cluster 3
        else
         con_r_c3 = [];   
        end
        
        if ~isempty(ripples_c4)
        idx7 = ismember(ripples(:,1:3),ripples_c4,'rows');
        con_r_c4 = ripples(idx7, 5:7);                   % Cluster 4
        else 
            con_r_c4 = [];
        end 
end 
%% Delta-spindle sequence detection
        seq_del_spin = [];
           if ~isnan(deltas) 
               if ~isnan(spindles)
                   co=[];
                   v_c1 = 0;
                   seq_indices =  find(co);   
                    win=0.1;
                    min_diff=100; %milliseconds
                    max_diff=1300; %milliseconds
                
                    for c = 1:length(deltas)
                        co = (spindles(:,7)>(deltas(c,2)+min_diff/1000)) & (spindles(:,7)<(deltas(c,2)+max_diff/1000));                   
                        if sum(co)~=0
                        seq_indices =  find(co);   
                              for f = 1 :length(seq_indices)
                                  v_c1 = v_c1+1;
                                  seq_del_spin(v_c1,1) = deltas(c,1);
                                  seq_del_spin(v_c1,2) = deltas(c,2);
                                  seq_del_spin(v_c1,3) = deltas(c,3);
                                  seq_del_spin(v_c1,4) = spindles(seq_indices(f),6);
                                  seq_del_spin(v_c1,5) = spindles(seq_indices(f),7);
                                  seq_del_spin(v_c1,6) = spindles(seq_indices(f),8);
                              end                         
                        end
                    end
               else 
                   seq_del_spin = NaN;
               end
           else 
               seq_del_spin = NaN;
           end
 %% Cooccurrence detection 
 
         s_start = spindles(:,6); % All Spindles
         s_peak = spindles(:,7);
         s_end = spindles(:,8);
         
         seq_s_start = seq_del_spin(:,4); % Sequential Spindles
         seq_s_peak = seq_del_spin(:,5);
         seq_s_end = seq_del_spin(:,6);
         
     if ~isempty(con_r_c1)
         r_start_c1 = con_r_c1(:,1);
         r_peak_c1 = con_r_c1(:,2);
         r_end_c1 = con_r_c1(:,3);
         
         [co_vec1, co_vec2, count_co_vec1, count_co_vec2] = cooccurrence_vec(r_start_c1,r_end_c1, s_start, s_end);% HPC, Cortex
         
         CO_ripples_c1{i} = co_vec1; % Indices 
         CO_spindles_c1{i} = co_vec2;

         CO_count_hpc_c1(i) = count_co_vec1; % Counts
         CO_count_pfc_c1(i) = count_co_vec2;
         
         [co_vec1, co_vec2, count_co_vec1, count_co_vec2] = cooccurrence_vec(r_start_c1,r_end_c1, seq_s_start, seq_s_end);% HPC, Cortex
 
         del_COspin_count_c1{i} = count_co_vec2;
         
     else 
         CO_ripples_c1{i} = NaN;
         CO_spindles_c1{i} = NaN;

         CO_count_hpc_c1(i) = 0;
         CO_count_pfc_c1(i) = 0;
         
         del_COspin_count_c1{i} = 0;
     end 
     if ~isempty(con_r_c2)
  
         r_start_c2 = con_r_c2(:,1);
         r_peak_c2 = con_r_c2(:,2);
         r_end_c2 = con_r_c2(:,3);
         
         [co_vec1, co_vec2, count_co_vec1, count_co_vec2] = cooccurrence_vec(r_start_c2,r_end_c2, s_start, s_end);% HPC, Cortex
         
         CO_ripples_c2{i} = co_vec1; % Indices 
         CO_spindles_c2{i} = co_vec2;

         CO_count_hpc_c2(i) = count_co_vec1; % Counts
         CO_count_pfc_c2(i) = count_co_vec2;
         
        [co_vec1, co_vec2, count_co_vec1, count_co_vec2] = cooccurrence_vec(r_start_c2,r_end_c2, seq_s_start, seq_s_end);% HPC, Cortex
 
         del_COspin_count_c2{i} = count_co_vec2;

     else 
         CO_ripples_c2{i} = NaN;
         CO_spindles_c2{i} = NaN;

         CO_count_hpc_c2(i) = 0;
         CO_count_pfc_c2(i) = 0;
         
         del_COspin_count_c2{i} = 0;
     end 
                  
     if ~isempty(con_r_c3)
         r_start_c3 = con_r_c3(:,1);
         r_peak_c3 = con_r_c3(:,2);
         r_end_c3 = con_r_c3(:,3);
         
         [co_vec1, co_vec2, count_co_vec1, count_co_vec2] = cooccurrence_vec(r_start_c3,r_end_c3, s_start, s_end);% HPC, Cortex
         
         CO_ripples_c3{i} = co_vec1; % Indices 
         CO_spindles_c3{i} = co_vec2;

         CO_count_hpc_c3(i) = count_co_vec1; % Counts
         CO_count_pfc_c3(i) = count_co_vec2;
         
        [co_vec1, co_vec2, count_co_vec1, count_co_vec2] = cooccurrence_vec(r_start_c3,r_end_c3, seq_s_start, seq_s_end);% HPC, Cortex
 
        del_COspin_count_c3{i} = count_co_vec2;

     else 
         CO_ripples_c3{i} = NaN;
         CO_spindles_c3{i} = NaN;

         CO_count_hpc_c3(i) = 0;
         CO_count_pfc_c3(i) = 0;
         
         del_COspin_count_c3{i} = 0;

     end 
     
     if ~isempty(con_r_c4)
         r_start_c4 = con_r_c4(:,1);
         r_peak_c4 = con_r_c4(:,2);
         r_end_c4 = con_r_c4(:,3);
         
         [co_vec1, co_vec2, count_co_vec1, count_co_vec2] = cooccurrence_vec(r_start_c4,r_end_c4, s_start, s_end);% HPC, Cortex
      
         CO_ripples_c4{i} = co_vec1; % Indices 
         CO_spindles_c4{i} = co_vec2;

         CO_count_hpc_c4(i) = count_co_vec1; % Counts
         CO_count_pfc_c4(i) = count_co_vec2;
         
         [co_vec1, co_vec2, count_co_vec1, count_co_vec2] = cooccurrence_vec(r_start_c4,r_end_c4, seq_s_start, seq_s_end);% HPC, Cortex
 
         del_COspin_count_c4{i} = count_co_vec2;
         
     else 
         CO_ripples_c4{i} = NaN;
         CO_spindles_c4{i} = NaN;

         CO_count_hpc_c4(i) = 0;
         CO_count_pfc_c4(i) = 0;
         
         del_COspin_count_c4{i} = 0;

     end 
         
         else
            
         CO_ripples_c1{i} = NaN;
         CO_spindles_c1{i} = NaN;

         CO_count_hpc_c1(i) = 0;
         CO_count_pfc_c1(i) = 0;
         
         CO_ripples_c2{i} = NaN;
         CO_spindles_c2{i} = NaN;

         CO_count_hpc_c2(i) = 0;
         CO_count_pfc_c2(i) = 0;
         
         CO_ripples_c3{i} = NaN;
         CO_spindles_c3{i} = NaN;

         CO_count_hpc_c3(i) = 0;
         CO_count_pfc_c3(i) = 0;
         
         CO_ripples_c4{i} = NaN;
         CO_spindles_c4{i} = NaN;

         CO_count_hpc_c4(i) = 0;
         CO_count_pfc_c4(i) = 0;
         
         del_COspin_count_c1{i} = 0; 
         del_COspin_count_c2{i} = 0;
         del_COspin_count_c3{i} = 0;
         del_COspin_count_c4{i} = 0;    
         
        end 
    end
    
save (strcat('Cluster1_CO_',SD_folders{j},'.mat'),'CO_ripples_c1','CO_spindles_c1','CO_count_hpc_c1','CO_count_pfc_c1')
save (strcat('Cluster2_CO_',SD_folders{j},'.mat'),'CO_ripples_c2','CO_spindles_c2','CO_count_hpc_c2','CO_count_pfc_c2')
save (strcat('Cluster3_CO_',SD_folders{j},'.mat'),'CO_ripples_c3','CO_spindles_c3','CO_count_hpc_c3','CO_count_pfc_c3')
save (strcat('Cluster4_CO_',SD_folders{j},'.mat'),'CO_ripples_c4','CO_spindles_c4','CO_count_hpc_c4','CO_count_pfc_c4')
save (strcat('Cluster1_4_Del_COspin_count_',SD_folders{j},'.mat'),'del_COspin_count_c1','del_COspin_count_c2','del_COspin_count_c3','del_COspin_count_c4')

end