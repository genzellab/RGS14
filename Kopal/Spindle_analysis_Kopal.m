clear variables
addpath(genpath('/Volumes/Samsung_T5/Milan_DA/OS_ephys_da/CorticoHippocampal'))
addpath('/Volumes/Samsung_T5/Milan_DA/RGS14_Ephys_da')
addpath ('/Volumes/Samsung_T5/Milan_DA/OS_ephys_da/ADRITOOLS')
addpath(genpath('/Volumes/Samsung_T5/Milan_DA/OS_ephys_da/FMAToolbox'))
cd('/Volumes/Samsung_T5/Milan_DA/RGS14_Ephys_da/Data_RGS14_Downsampled_First_Session')
% cd('/home/adrian/Documents/new_OS_downsampled_2020')

% ----spindles_bout_specific_timestamps ---> Misleading variable name;
% These are the ones wrt to the signal and have the same format as the
% ripple ones; they don't restart per nrem bout.

%Select sleep stage or use all.

list = {'Wake','NREM','Transitional','REM','All'};
[indx1] = listdlg('SelectionMode','single','ListString',list,'InitialValue',2);

switch indx1
    case 1
        ss = 1; 
        stage='Wake';
    case 2
        ss = 3; 
        stage='NREM';
    case 3
        ss = 4; %Transitional
        stage='Trans';
    case 4
        ss = 5; %REM
        stage='REM';
    case 5
        ss = 6; %ALL
        stage='ALL';
        
end

   
fn=1000;

total_spindles=[];
total_spindles_minute=[];
phase_total_complete={};
rat_folder=getfolder;
% for k=1:length(rat_folder)
prompt = {'Enter the rat index'};
    dlgtitle = 'Rat Index';
    % opts.Interpreter = 'tex';
  k = str2double(inputdlg(prompt,dlgtitle));
% for k=1 % rat index 
    cd(rat_folder{k})    
    g=getfolder;
    spindle_phase_comp = [];
    spindles_waveforms_comp = [];
    spindles_waveforms_broadband_comp = [];

for j=1:length(g) %study day index
%     for j=3

       cd(g{j})
%end
        G=getfolder;
                
%%        
%Get presleep
cfold3=[];
cfold=G(or(cellfun(@(x) ~isempty(strfind(x,'pre')),G),cellfun(@(x) ~isempty(strfind(x,'Pre')),G)));
for q=1:length(cfold)
    if (~contains(cfold{q}, 'test') && ~contains(cfold{q}, 'Test'))
        cfold3=[cfold3; cfold{q}];
    end
end
if ~isempty(cfold3)
    cfold=cellstr(cfold3)';
end

% Get post trials
cfold3=[];
cfold2=G(or(cellfun(@(x) ~isempty(strfind(x,'post')),G),cellfun(@(x) ~isempty(strfind(x,'Post')),G)));
for q=1:length(cfold2)
    if (~contains(cfold2{q}, 'test') && ~contains(cfold2{q}, 'Test'))
        cfold3=[cfold3; cfold2{q}];
    end
end
cfold2=cellstr(cfold3)';

%Get presleep
% cfold=G(or(cellfun(@(x) ~isempty(strfind(x,'pre')),G),cellfun(@(x) ~isempty(strfind(x,'Pre')),G)));
% cfold=cfold(and(~contains(cfold,'test'),~contains(cfold,'Test')));
% 
% % Get post trials
% cfold2=G(or(cellfun(@(x) ~isempty(strfind(x,'post')),G),cellfun(@(x) ~isempty(strfind(x,'Post')),G)));
% cfold2=cfold2(and(~contains(cfold2,'test'),~contains(cfold2,'Test')));
%%
%Ignore trial 6
for ind=1:length(cfold2)
  if  ~(contains(cfold2{ind},'trial1') ||contains(cfold2{ind},'trial2')||contains(cfold2{ind},'trial3')||contains(cfold2{ind},'trial4')||contains(cfold2{ind},'trial5')...
        ||contains(cfold2{ind},'Trial1')||contains(cfold2{ind},'Trial2')||contains(cfold2{ind},'Trial3')||contains(cfold2{ind},'Trial4')||contains(cfold2{ind},'Trial5')  )
      
      cfold2{ind}=[];    
  end
end

cfold2=cfold2(~cellfun('isempty',cfold2));

G=[cfold cfold2];

        
              
     if isempty(G)
            no_folder=1;
            %g=NaN;
        else
            no_folder=0;

         for i=1:length(G); %trial index 
%          for i=6
                clear states
%             for i=1
%                  xo
                cd(G{i})
                A = dir('*states*.mat');
                A={A.name};
                
                if sum(contains(A, 'states')) > 0 %More than 2 sleep scoring files
                    A=A(cellfun(@(x) ~isempty(strfind(x,'states')),A));
                    A=A(~(cellfun(@(x) ~isempty(strfind(x,'eeg')),A)));
                    
                    if sum(contains(A, 'states')) > 0
%                     st2=st(cellfun(@(x) ~isempty(strfind(x,barea)),st)); %Brain area.
                    cellfun(@load,A);



Cortex=dir(strcat('*','PFC','*.mat'));
Cortex=Cortex.name;
Cortex=load(Cortex);
% Cortex=Cortex.Cortex;
Cortex=getfield(Cortex,'PFC');
Cortex = Cortex.*(0.195);


% Wn1=[.5/(fn/2) 4/(fn/2)]; % Cutoff=.5-4 Hz
% % Wn1=[50/(fn/2) 80/(fn/2)]; 
% [b3,a3] = butter(3,Wn1,'bandpass'); %Filter coefficients
% Cortex_delta= filtfilt(b3,a3,Cortex);
% cortex_phase=mod(rad2deg(angle(hilbert(Cortex_delta))),360);

% continuous_timestamps=[(0:length(Cortex)-1)/1000];


                                      

                    if and(~contains(G{i},'trial5'),~contains(G{i},'Trial5')) %Whenever it is not PostTrial 5 
                        
                        % Sleep Scoring data
                        if length(states)<45*60
                            states=[states nan(1,45*60-length(states))]; %Fill with NaNs.
                        else
                            states=states(1:45*60); %Take only 45 min.
                        end
                        
                        %Ephys data
                        
                        if length(Cortex)<45*60*1000
                            Cortex=[Cortex.' (nan(45*60*1000-length(Cortex),1).')]; %Fill with NaNs.
                        else
                            Cortex=Cortex(1:45*60*1000).'; %Take only 45 min.
                        end
                         
                    continuous_timestamps=[(0:length(Cortex)-1)/1000];
                        %Ignore NaNs
                        PFC=Cortex;
    if sum(isnan(PFC))~=0 
        PFC(isnan(PFC)) = 0;
        states(isnan(states)) = 0;
    end



%LPF 300 Hz:
Wn1=[320/(fn/2)]; % Cutoff=320 Hz
[b2,a2] = butter(3,Wn1); %Filter coefficients


%Convert signal to 1 sec epochs.
    e_t=1;
    e_samples=e_t*(1000); %fs=1kHz
    ch=length(PFC);
    nc=floor(ch/e_samples); %Number of epochs
    NC=[];
    NC2=[];
    NC3=[];
    
    for kk=1:nc    

      NC2(:,kk) = PFC(1+e_samples*(kk-1):e_samples*kk);
      NC3(:,kk) = continuous_timestamps(1+e_samples*(kk-1):e_samples*kk);
    end
    
    vec_bin=states;
    vec_bin(vec_bin~=ss)=0;
    vec_bin(vec_bin==ss)=1;
    
    
  %% Event Detection 
   [spindles,phase_total,z,BST,V_pfc,v_pfc]= spindles_detection(vec_bin,NC3,NC2,fn,j,i,k);
   
   if iscell(V_pfc) && ~isempty(spindles)
    concatenated_NREM_pfc = vertcat(V_pfc{:});
    waveforms_spindles={};
        for c=1:size(spindles,1)
           waveforms_spindles{c,1}= concatenated_NREM_pfc(int32(spindles(c,6)*1000+1):int32(spindles(c,8)*1000+1));
        end
    else
        waveforms_spindles = NaN;
   end 
   
   if iscell(v_pfc) && ~isempty(spindles)
    concatenated_NREM_broadband_pfc = vertcat(v_pfc{:});
    waveforms_spindles_broadband={};
        for c=1:size(spindles,1)
           waveforms_spindles_broadband{c,1}= concatenated_NREM_broadband_pfc(int32(spindles(c,6)*1000+1):int32(spindles(c,8)*1000+1));
        end
    else
        waveforms_spindles_broadband = NaN;
   end 
   spindles_waveform_broadband_total{i} = waveforms_spindles_broadband;
   spindles_waveform_total{i} = waveforms_spindles;
   spindles_complete{j,i} = spindles;
   spindles_bout_specific_timestamps{i} = BST;
   phase_total_complete{j,1}{i,1} = phase_total;
   total_spindles(j,i) = z;
   stage_count = sum(states(:)==ss);
   total_spindles_minute = (z/stage_count*60);
   total_spindles_minute_complete(j,i) = total_spindles_minute;
    
clear V_pfc v_pfc V_pfc_bp V_pfc_bp2 v_values vec_bin VV_pfc v2 v_index NC2 NC3 VV_timestamps VV_pfc_phase v_timestamps
                    
                                  
                    elseif contains(G{i}, 'rial5') % Post Trial 5 case 
% %                         xo

                        %Sleep scoring data
                        
                        if length(states)<45*60*4
                            states = [states nan(1,45*60*4-length(states))]; %Fill with NaNs.
                        else
                            states = states(1:45*60*4); %Take only 45 min.
                        end
                        
                        
                        %Ephys

                        
                        if length(Cortex)<45*60*1000*4
                            Cortex = [Cortex.' (nan(45*60*1000*4-length(Cortex),1).')]; %Fill with NaNs.
                        else
                            Cortex = Cortex(1:45*60*1000*4).'; %Take only 45 min.
                        end

                        for jj=1:4
% for jj=3
                            
                                PFC = Cortex((2700*1000*(jj-1))+1:2700*1000*jj);
                                states_chunk = states(2700*(jj-1)+1:2700*jj);
                                continuous_timestamps = [(0:length(PFC)-1)/1000];
                                if sum(isnan(PFC))~=0 
                                PFC(isnan(PFC)) = 0;
                                states_chunk(isnan(states_chunk)) = 0;
                                end
                                    e_t = 1;
                                    e_samples = e_t*(1000); %fs=1kHz
                                    ch = length(PFC);
                                    nc = floor(ch/e_samples); %Number of epochs  
                                    NC2 = [];
                                    NC3 = [];
                                        for kk=1:nc    

                                          NC2(:,kk) = PFC(1+e_samples*(kk-1):e_samples*kk);
                                          NC3(:,kk) = continuous_timestamps(1+e_samples*(kk-1):e_samples*kk);

                                        end

                                    vec_bin = states_chunk;
                                    vec_bin(vec_bin~=ss)=0;
                                    vec_bin(vec_bin==ss)=1; 
%% Event Detection                                     
                                   [spindles,phase_total,z,BST,V_pfc,v_pfc] = spindles_detection(vec_bin,NC3,NC2,fn,j,i,k);
                                   if iscell(V_pfc) && ~isempty(spindles)
                                    concatenated_NREM_pfc = vertcat(V_pfc{:});
                                    waveforms_spindles = {};
                                        for c = 1:size(spindles,1)
                                           waveforms_spindles{c,1} = concatenated_NREM_pfc(int32(spindles(c,6)*1000+1):int32(spindles(c,8)*1000+1));
                                        end
                                    else
                                        waveforms_spindles = NaN;
                                   end 
                                   if iscell(v_pfc) && ~isempty(spindles)
                                    concatenated_NREM_broadband_pfc = vertcat(v_pfc{:});
                                    waveforms_spindles_broadband = {};
                                        for c = 1:size(spindles,1)
                                           waveforms_spindles_broadband{c,1} = concatenated_NREM_broadband_pfc(int32(spindles(c,6)*1000+1):int32(spindles(c,8)*1000+1));
                                        end
                                    else
                                        waveforms_spindles_broadband = NaN;
                                   end 
                                   spindles_waveform_broadband_total{i+jj-1}=waveforms_spindles_broadband;
                                   spindles_waveform_total{i+jj-1}=waveforms_spindles;
                                   spindles_complete{j,i+jj-1}= spindles;
                                   spindles_bout_specific_timestamps{i+jj-1}= BST;
                                   phase_total_complete{j,1}{i+jj-1,1} = phase_total;
                                   total_spindles(j,i+jj-1) = z;
                                   stage_count = sum(states(:)==ss);
                                   total_spindles_minute = (z/stage_count*60);
                                   total_spindles_minute_complete(j,i+jj-1) = total_spindles_minute;
clear V_pfc v_pfc V_pfc_bp V_pfc_bp2 v_values vec_bin VV_pfc v2 v_index NC2 NC3
                            
                       end 
                    else
                        continue
                    end
                else
                    cd ..
                    continue
                end
                cd ..
            end
        end
        cd ..
% %         spindles_total_data = spindles_total_data_complete{j};

        phase_total = phase_total_complete{j}.';
        spindle_phase_comp = [spindle_phase_comp;phase_total];
        spindles_count = total_spindles(j,:);
        spindles_waveforms_comp = [spindles_waveforms_comp; spindles_waveform_total];
        spindles_waveforms_broadband_comp = [spindles_waveforms_broadband_comp; spindles_waveform_broadband_total];

        spindles_total_data = spindles_complete(j,:);
        
        save(strcat('spindles_waveforms_',g{j},'.mat'),'spindles_waveform_total','-v7.3')
        save(strcat('spindles_waveforms_broadband_',g{j},'.mat'),'spindles_waveform_broadband_total','-v7.3')
        
        save(strcat('spindles_count_',g{j},'.mat'),'spindles_count')

        save(strcat('spindles_phase_',g{j},'.mat'),'phase_total')
        save(strcat('spindles_timestamps_',g{j},'.mat'),'spindles_bout_specific_timestamps') 
        save(strcat('spindles_total_data_',g{j},'.mat'),'spindles_total_data') 
        
        if j == length(g)
        save(strcat('spindles_data_compilation_Rat',rat_folder{k},'.mat'),'total_spindles','total_spindles_minute','spindles_complete')
        end
        
        save(strcat('spindles_phases_compilation_Rat',rat_folder{k},'.mat'),'spindle_phase_comp' )
        save(strcat('spindles_waveforms_compilation_Rat',rat_folder{k},'.mat'),'spindles_waveforms_comp','-v7.3')
        save(strcat('spindles_waveforms_broadband_compilation_Rat',rat_folder{k},'.mat'),'spindles_waveforms_broadband_comp','-v7.3')

    end 
end 

