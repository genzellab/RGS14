clear variables
addpath '/Volumes/Samsung_T5/Milan_DA/RGS14_Ephys_data'
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
CO_ripples = {};
CO_spindles = {};
CO_count_hpc = [];
CO_count_pfc = [];

%% Extraction of indices and counts wrt Co-occurring events
    for i = 1:length(spindles_total_data)
        if ~isnan(spindles_total_data{1,i})
         s_start = spindles_total_data{1,i}(:,6);
         s_peak = spindles_total_data{1,i}(:,7);
         s_end = spindles_total_data{1,i}(:,8);

         r_start = ripple_total_data{1,i}(:,5);
         r_peak = ripple_total_data{1,i}(:,6);
         r_end = ripple_total_data{1,i}(:,7);
         
 [co_vec1, co_vec2, count_co_vec1, count_co_vec2] = cooccurrence_vec(r_start,r_end, s_start, s_end);% HPC, Cortex

         CO_ripples{i} = co_vec1; % Indices 
         CO_spindles{i} = co_vec2;

         CO_count_hpc(i) = count_co_vec1; % Counts
         CO_count_pfc(i) = count_co_vec2;

        else
            
         CO_ripples{i} = NaN;
         CO_spindles{i} = NaN;

         CO_count_hpc(i) = NaN;
         CO_count_pfc(i) = NaN;
        end 
    end
    save (strcat('CO_',SD_folders{j},'.mat'),'CO_ripples','CO_spindles','CO_count_hpc','CO_count_pfc')
end