clear variables

addpath(genpath('/Volumes/Samsung_T5/Milan_DA/OS_ephys_da/CorticoHippocampal'))
addpath('/Volumes/Samsung_T5/Milan_DA/RGS14_Ephys_da')
addpath ('/Volumes/Samsung_T5/Milan_DA/OS_ephys_da/ADRITOOLS')
addpath(genpath('/Volumes/Samsung_T5/Milan_DA/OS_ephys_da/FMAToolbox'))
cd('/Volumes/Samsung_T5/Milan_DA/RGS14_Ephys_da/Data_RGS14_Downsampled_First_Session')

% cd('/media/genzel/Data/Milan_DA/RGS14_Ephys_da/Data_RGS14_Downsampled')
% cd('/home/adrian/Documents/new_OS_downsampled_2020')

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

total_delta=[];
total_delta_minute=[];
delta_total_data={};
rat_folder=getfolder;
% for k=1:length(rat_folder)
for k = 1
    cd(rat_folder{k})    
    g=getfolder;
delta_waveform_broadband_comp = [];
delta_waveform_comp = [];
    for j=1:length(g)
    delta_waveform_broadband_total = {};
    delta_waveform_total = {};
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

            for i=1:length(G)
            
%          for i=1
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
Cortex=Cortex.*(0.195);



                                      

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
    
    for kk=1:nc    

      NC2(:,kk)= PFC(1+e_samples*(kk-1):e_samples*kk);
    end
    
    vec_bin=states;
    vec_bin(vec_bin~=ss)=0;
    vec_bin(vec_bin==ss)=1;
    
    if sum(vec_bin)==0  %%All states
          delta_total_data{j,1}{i,1} = NaN;
        if ss==6  %Only continue if All states option was selected
        vec_bin=vec_bin+1;
        else  %When the sleep state was not found.


V_pfc=0;
         cd ..
         continue
        end
    end
    %Cluster one values:
    v2=ConsecutiveOnes(vec_bin);
    
    v_index=find(v2~=0);
    v_values=v2(v2~=0);

%     
%     ver=NC(:, v_index(1):v_index(1)+(v_values(1,1)-1));
%     v{1}=reshape(A, numel(A), 1);

for epoch_count = 1:length(v_index)
v_pfc{epoch_count,1}=reshape(NC2(:, v_index(epoch_count):v_index(epoch_count)+(v_values(1,epoch_count)-1)), [], 1);
end 

%v_hpc and v_pfc: NREM epochs.

%Ripple detection

Wn1=[1/(fn/2) 6/(fn/2)]; % Cutoff=320 Hz
[b2,a2] = butter(3,Wn1); %Filter coefficients
V_pfc = cellfun(@(equis)filtfilt(b2,a2,equis), v_pfc ,'UniformOutput',false);
%LPF 300 Hz:
VV_pfc=[];
for b=1:size(V_pfc,1)
    VV_pfc=[VV_pfc; V_pfc{b}];
end
V_pfc=VV_pfc;
V_pfc={V_pfc};
z=0;
for p=1:length(V_pfc)
    V_pfc_bp = V_pfc{p};
    V_pfc_bp2=[];
    for l=1:length(V_pfc_bp)
        V_pfc_bp2(l,1)=l/1000;
    end
    V_pfc_bp=horzcat(V_pfc_bp2, V_pfc_bp);
    if length(V_pfc_bp2)>4000
        delta = FindDeltaWavesRGS14(V_pfc_bp,k);
        delta_total_data{j,1}{i,p} = delta;
        z=z+size(delta,1); 
     else
     delta_total_data{j,1}{i,p} = NaN;
    end
    
end

if iscell(v_pfc) && ~isempty(delta)
    concatenated_NREM_broadband_pfc = vertcat(v_pfc{:});
    waveforms_delta_broadband={};
        for c=1:size(delta,1)
           waveforms_delta_broadband{c,1}= concatenated_NREM_broadband_pfc(int32(delta(c,1)*1000+1):int32(delta(c,3)*1000+1));
        end
else
        waveforms_delta_broadband = NaN;
end
if iscell(V_pfc) && ~isempty(delta)
    concatenated_NREM_pfc = vertcat(V_pfc{:});
    waveforms_delta={};
        for c=1:size(delta,1)
           waveforms_delta{c,1}= concatenated_NREM_pfc(int32(delta(c,1)*1000+1):int32(delta(c,3)*1000+1));
        end
else
        waveforms_delta = NaN;
end
delta_waveform_broadband_total{i} = waveforms_delta_broadband;
delta_waveform_total{i} = waveforms_delta;  
total_delta(j,i)=z;
stage_count = sum(states(:)==ss);
total_delta_minute(j,i)=(total_delta(j,i)/stage_count*60);

clear V_pfc v_pfc V_pfc_bp V_pfc_bp2 v_values vec_bin VV_pfc v2 v_index NC2 NC

                    elseif contains(G{i}, 'rial5') % PostTrial 5 case 
% %                         xo
                        %Sleep scoring data
                        if length(states)<45*60*4
                            states=[states nan(1,45*60*4-length(states))]; %Fill with NaNs.
                        else
                            states=states(1:45*60*4); %Take only 45 min.
                        end
                        
                        
                        %Ephys

                        
                        if length(Cortex)<45*60*1000*4
                            Cortex=[Cortex.' (nan(45*60*1000*4-length(Cortex),1).')]; %Fill with NaNs.
                        else
                            Cortex=Cortex(1:45*60*1000*4).'; %Take only 45 min.
                        end

                        
                        %Chunk in 4
                        states1=states(1:2700);
                        states2=states(2700+1:2700*2);
                        states3=states(1+2700*2:2700*3);
                        states4=states(1+2700*3:2700*4);
                    
                        
                        PFC_1=Cortex(1:2700*1000);
                        PFC_2=Cortex(2700*1000+1:2700*2*1000);
                        PFC_3=Cortex(1+2700*2*1000:2700*3*1000);
                        PFC_4=Cortex(1+2700*3*1000:2700*4*1000);
                        
                            if sum(isnan(PFC_1))~=0 
                                PFC_1(isnan(PFC_1)) = 0;
                                states1(isnan(states1)) = 0;
                            end
                                                        if sum(isnan(PFC_2))~=0 
                                PFC_2(isnan(PFC_2)) = 0;
                                states2(isnan(states2)) = 0;
                                                        end
                                                        if sum(isnan(PFC_3))~=0 
                                PFC_3(isnan(PFC_3)) = 0;
                                states3(isnan(states3)) = 0;
                                                        end
                                                        if sum(isnan(PFC_4))~=0 
                                PFC_4(isnan(PFC_4)) = 0;
                                states4(isnan(states4)) = 0;
                            end
                        
                            %LPF 300 Hz:
Wn1=[320/(fn/2)]; % Cutoff=320 Hz
[b2,a2] = butter(3,Wn1); %Filter coefficients


    e_t=1;
    e_samples=e_t*(1000); %fs=1kHz
    ch=length(PFC_1);
    nc=floor(ch/e_samples); %Number of epochs
    NC=[];
    NC2=[];
    
    for kk=1:nc    

      NC2(:,kk)= PFC_1(1+e_samples*(kk-1):e_samples*kk);
    end
    
    vec_bin=states1;
    vec_bin(vec_bin~=ss)=0;
    vec_bin(vec_bin==ss)=1;
    
    if sum(vec_bin)~=0  %%All states
%         if ss==6  %Only continue if All states option was selected
%         vec_bin=vec_bin+1;
%         else  %When the sleep state was not found.
% 
% 
% V_pfc=0;
%          cd ..
%          continue
%         end
%     end
    %Cluster one values:
    v2=ConsecutiveOnes(vec_bin);
    
    v_index=find(v2~=0);
    v_values=v2(v2~=0);

%     
%     ver=NC(:, v_index(1):v_index(1)+(v_values(1,1)-1));
%     v{1}=reshape(A, numel(A), 1);
for epoch_count=1:length(v_index)

v_pfc{epoch_count,1}=reshape(NC2(:, v_index(epoch_count):v_index(epoch_count)+(v_values(1,epoch_count)-1)), [], 1);
end 

%v_hpc and v_pfc: NREM epochs.

%Ripple detection

Wn1=[1/(fn/2) 6/(fn/2)]; % Cutoff=320 Hz
[b2,a2] = butter(3,Wn1); %Filter coefficients
V_pfc=cellfun(@(equis) filtfilt(b2,a2,equis), v_pfc ,'UniformOutput',false);
%LPF 300 Hz:
VV_pfc=[];
for b=1:size(V_pfc,1)
    VV_pfc=[VV_pfc; V_pfc{b}];
end
V_pfc=VV_pfc;
V_pfc={V_pfc};
z=0;
for p=1:length(V_pfc)
    V_pfc_bp = V_pfc{p};
    V_pfc_bp2=[];
    for l=1:length(V_pfc_bp)
        V_pfc_bp2(l,1)=l/1000;
    end
    V_pfc_bp=horzcat(V_pfc_bp2, V_pfc_bp);
    if length(V_pfc_bp2)>4000
        delta = FindDeltaWavesRGS14(V_pfc_bp,k);
        delta_total_data{j,1}{i,p} = delta;
        z=z+size(delta,1); 
     else
     delta_total_data{j,1}{i,p} = NaN;
    end
    
end
if iscell(v_pfc) && ~isempty(delta)
    concatenated_NREM_broadband_pfc = vertcat(v_pfc{:});
    waveforms_delta_broadband={};
        for c=1:size(delta,1)
           waveforms_delta_broadband{c,1}= concatenated_NREM_broadband_pfc(int32(delta(c,1)*1000+1):int32(delta(c,3)*1000+1));
        end
    else
        waveforms_delta_broadband = NaN;
end
if iscell(V_pfc) && ~isempty(delta)
    concatenated_NREM_pfc = vertcat(V_pfc{:});
    waveforms_delta={};
        for c=1:size(delta,1)
           waveforms_delta{c,1}= concatenated_NREM_pfc(int32(delta(c,1)*1000+1):int32(delta(c,3)*1000+1));
        end
else
        waveforms_delta = NaN;
end
delta_waveform_broadband_total{i} = waveforms_delta_broadband; 
delta_waveform_total{i} = waveforms_delta;
total_delta(j,i) = z;
stage_count = sum(states1(:)==ss);
total_delta_minute(j,i) = (total_delta(j,i)/stage_count*60);
else 
      delta_total_data{j,1}{i,p} = NaN;
      delta_waveform_broadband_total{i} = [];
      delta_waveform_total{i} = [];
   end
    
clear V_pfc v_pfc V_pfc_bp V_pfc_bp2 v_values vec_bin VV_pfc v2 v_index NC2 NC

    e_t=1;
    e_samples=e_t*(1000); %fs=1kHz
    ch=length(PFC_2);
    nc=floor(ch/e_samples); %Number of epochs
    NC=[];
    NC2=[];
    
    for kk=1:nc    

      NC2(:,kk)= PFC_2(1+e_samples*(kk-1):e_samples*kk);
    end
    
    vec_bin=states2;
    vec_bin(vec_bin~=ss)=0;
    vec_bin(vec_bin==ss)=1;
    
    if sum(vec_bin)~=0  %%All states
%         if ss==6  %Only continue if All states option was selected
%         vec_bin=vec_bin+1;
%         else  %When the sleep state was not found.
% 
% 
% V_pfc=0;
%          cd ..
%          continue
%         end
%     end
    %Cluster one values:
    v2=ConsecutiveOnes(vec_bin);
    
    v_index=find(v2~=0);
    v_values=v2(v2~=0);

%     
%     ver=NC(:, v_index(1):v_index(1)+(v_values(1,1)-1));
%     v{1}=reshape(A, numel(A), 1);
for epoch_count=1:length(v_index)

v_pfc{epoch_count,1}=reshape(NC2(:, v_index(epoch_count):v_index(epoch_count)+(v_values(1,epoch_count)-1)), [], 1);

end 

%v_hpc and v_pfc: NREM epochs.

%Ripple detection

Wn1=[1/(fn/2) 6/(fn/2)]; % Cutoff=320 Hz
[b2,a2] = butter(3,Wn1); %Filter coefficients
V_pfc=cellfun(@(equis) filtfilt(b2,a2,equis), v_pfc ,'UniformOutput',false);
%LPF 300 Hz:
VV_pfc=[];
for b=1:size(V_pfc,1)
    VV_pfc=[VV_pfc; V_pfc{b}];
end
V_pfc=VV_pfc;
V_pfc={V_pfc};
z=0;
for p=1:length(V_pfc)
    V_pfc_bp = V_pfc{p};
    V_pfc_bp2=[];
    for l=1:length(V_pfc_bp)
        V_pfc_bp2(l,1)=l/1000;
    end
    V_pfc_bp=horzcat(V_pfc_bp2, V_pfc_bp);
    if length(V_pfc_bp2)>4000
        delta = FindDeltaWavesRGS14(V_pfc_bp,k);
        delta_total_data{j,1}{i+1,p} = delta;
        z=z+size(delta,1);      
   else
         delta_total_data{j,1}{i+1,p} = NaN;
    end
    
end
if iscell(v_pfc) && ~isempty(delta)
    concatenated_NREM_broadband_pfc = vertcat(v_pfc{:});
    waveforms_delta_broadband={};
        for c=1:size(delta,1)
           waveforms_delta_broadband{c,1}= concatenated_NREM_broadband_pfc(int32(delta(c,1)*1000+1):int32(delta(c,3)*1000+1));
        end
    else
        waveforms_delta_broadband = NaN;
end
if iscell(V_pfc) && ~isempty(delta)
    concatenated_NREM_pfc = vertcat(V_pfc{:});
    waveforms_delta={};
        for c=1:size(delta,1)
           waveforms_delta{c,1}= concatenated_NREM_pfc(int32(delta(c,1)*1000+1):int32(delta(c,3)*1000+1));
        end
else
        waveforms_delta = NaN;
end

delta_waveform_broadband_total{i+1} = waveforms_delta_broadband;
delta_waveform_total{i+1}  = waveforms_delta;
total_delta(j,i+1)=z;
stage_count = sum(states2(:)==ss);
total_delta_minute(j,i+1)=(total_delta(j,i+1)/stage_count*60);
    else 
      delta_total_data{j,1}{i+1,p} = NaN;
      delta_waveform_broadband_total{i+1} = [];
      delta_waveform_total{i+1} = [];
    end
clear V_pfc v_pfc V_pfc_bp V_pfc_bp2 v_values vec_bin VV_pfc v2 v_index NC2 NC

    e_t=1;
    e_samples=e_t*(1000); %fs=1kHz
    ch=length(PFC_3);
    nc=floor(ch/e_samples); %Number of epochs
    NC=[];
    NC2=[];
    
    for kk=1:nc    

      NC2(:,kk)= PFC_3(1+e_samples*(kk-1):e_samples*kk);
    end
    
    vec_bin=states3;
    vec_bin(vec_bin~=ss)=0;
    vec_bin(vec_bin==ss)=1;
    
    if sum(vec_bin)~=0  %%All states
%         if ss==6  %Only continue if All states option was selected
%         vec_bin=vec_bin+1;
%         else  %When the sleep state was not found.
% 
% 
% V_pfc=0;
%          cd ..
%          continue
%         end
%     end
    %Cluster one values:
    v2=ConsecutiveOnes(vec_bin);
    
    v_index=find(v2~=0);
    v_values=v2(v2~=0);

%     
%     ver=NC(:, v_index(1):v_index(1)+(v_values(1,1)-1));
%     v{1}=reshape(A, numel(A), 1);
for epoch_count=1:length(v_index)

v_pfc{epoch_count,1}=reshape(NC2(:, v_index(epoch_count):v_index(epoch_count)+(v_values(1,epoch_count)-1)), [], 1);
end 

%v_hpc and v_pfc: NREM epochs.

%Ripple detection

Wn1=[1/(fn/2) 6/(fn/2)]; % Cutoff=320 Hz
[b2,a2] = butter(3,Wn1); %Filter coefficients
V_pfc=cellfun(@(equis) filtfilt(b2,a2,equis), v_pfc ,'UniformOutput',false);
%LPF 300 Hz:
VV_pfc=[];
for b=1:size(V_pfc,1)
    VV_pfc=[VV_pfc; V_pfc{b}];
end
V_pfc=VV_pfc;
V_pfc={V_pfc};
z=0;
for p=1:length(V_pfc)
    V_pfc_bp = V_pfc{p};
    V_pfc_bp2=[];
    for l=1:length(V_pfc_bp)
        V_pfc_bp2(l,1)=l/1000;
    end
    V_pfc_bp=horzcat(V_pfc_bp2, V_pfc_bp);
    if length(V_pfc_bp2)>4000
        delta = FindDeltaWavesRGS14(V_pfc_bp,k);
        delta_total_data{j,1}{i+2,p} = delta;
        z=z+size(delta,1);
    else
         delta_total_data{j,1}{i+2,p} = NaN;
    end
    
end
if iscell(v_pfc) && ~isempty(delta)
    concatenated_NREM_broadband_pfc = vertcat(v_pfc{:});
    waveforms_delta_broadband={};
        for c=1:size(delta,1)
           waveforms_delta_broadband{c,1}= concatenated_NREM_broadband_pfc(int32(delta(c,1)*1000+1):int32(delta(c,3)*1000+1));
        end
    else
        waveforms_delta_broadband = NaN;
end
if iscell(V_pfc) && ~isempty(delta)
    concatenated_NREM_pfc = vertcat(V_pfc{:});
    waveforms_delta={};
        for c=1:size(delta,1)
           waveforms_delta{c,1}= concatenated_NREM_pfc(int32(delta(c,1)*1000+1):int32(delta(c,3)*1000+1));
        end
else
        waveforms_delta = NaN;
end
delta_waveform_broadband_total{i+2} = waveforms_delta_broadband;
delta_waveform_total{i+2} = waveforms_delta;
total_delta(j,i+2)=z;
stage_count = sum(states3(:)==ss);
total_delta_minute(j,i+2)=(total_delta(j,i+2)/stage_count*60);
    else 
      delta_total_data{j,1}{i+2,p} = NaN;
      delta_waveform_broadband_total{i+2} =[];
      delta_waveform_total{i+2} = [];
    end
clear V_pfc v_pfc V_pfc_bp V_pfc_bp2 v_values vec_bin VV_pfc v2 v_index NC2 NC

    e_t=1;
    e_samples=e_t*(1000); %fs=1kHz
    ch=length(PFC_4);
    nc=floor(ch/e_samples); %Number of epochs
    NC=[];
    NC2=[];
    
    for kk=1:nc    

      NC2(:,kk)= PFC_4(1+e_samples*(kk-1):e_samples*kk);
    end
    
    vec_bin=states4;
    vec_bin(vec_bin~=ss)=0;
    vec_bin(vec_bin==ss)=1;
    
    if sum(vec_bin)~=0  %%All states
%         if ss==6  %Only continue if All states option was selected
%         vec_bin=vec_bin+1;
%         else  %When the sleep state was not found.
% 
% 
% V_pfc=0;
%          cd ..
%          continue
%         end
%     end
    %Cluster one values:
    v2=ConsecutiveOnes(vec_bin);
    
    v_index=find(v2~=0);
    v_values=v2(v2~=0);

%     
%     ver=NC(:, v_index(1):v_index(1)+(v_values(1,1)-1));
%     v{1}=reshape(A, numel(A), 1);
for epoch_count=1:length(v_index)

v_pfc{epoch_count,1}=reshape(NC2(:, v_index(epoch_count):v_index(epoch_count)+(v_values(1,epoch_count)-1)), [], 1);

end 

%v_hpc and v_pfc: NREM epochs.

%Ripple detection

Wn1=[1/(fn/2) 6/(fn/2)]; % Cutoff=320 Hz
[b2,a2] = butter(3,Wn1); %Filter coefficients
V_pfc=cellfun(@(equis) filtfilt(b2,a2,equis), v_pfc ,'UniformOutput',false);
%LPF 300 Hz:
VV_pfc=[];
for b=1:size(V_pfc,1)
    VV_pfc=[VV_pfc; V_pfc{b}];
end
V_pfc=VV_pfc;
V_pfc={V_pfc};
z=0;
for p=1:length(V_pfc)
    V_pfc_bp = V_pfc{p};
    V_pfc_bp2=[];
    for l=1:length(V_pfc_bp)
        V_pfc_bp2(l,1)=l/1000;
    end
    V_pfc_bp=horzcat(V_pfc_bp2, V_pfc_bp);
    if length(V_pfc_bp2)>4000
        delta = FindDeltaWavesRGS14(V_pfc_bp,k);
        delta_total_data{j,1}{i+3,p} = delta;
        z=z+size(delta,1); 
    else
         delta_total_data{j,1}{i+3,p} = NaN;     
    end
    
end
if iscell(v_pfc) && ~isempty(delta)
    concatenated_NREM_broadband_pfc = vertcat(v_pfc{:});
    waveforms_delta_broadband={};
        for c=1:size(delta,1)
           waveforms_delta_broadband{c,1}= concatenated_NREM_broadband_pfc(int32(delta(c,1)*1000+1):int32(delta(c,3)*1000+1));
        end
    else
        waveforms_delta_broadband = NaN;
end
if iscell(V_pfc) && ~isempty(delta)
    concatenated_NREM_pfc = vertcat(V_pfc{:});
    waveforms_delta={};
        for c=1:size(delta,1)
           waveforms_delta{c,1}= concatenated_NREM_pfc(int32(delta(c,1)*1000+1):int32(delta(c,3)*1000+1));
        end
else
        waveforms_delta = NaN;
end
delta_waveform_broadband_total{i+3} = waveforms_delta_broadband; 
delta_waveform_total{i+3} = waveforms_delta; 
total_delta(j,i+3)=z;
stage_count = sum(states4(:)==ss);
total_delta_minute(j,i+3)=(total_delta(j,i+3)/stage_count*60);
    else 
      delta_total_data{j,1}{i+3,p} = NaN;
      delta_waveform_broadband_total{i+3} = [];
      delta_waveform_total{i+3} =[];
    end
clear V_pfc v_pfc V_pfc_bp V_pfc_bp2 v_values vec_bin VV_pfc v2 v_index NC2 NC
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
  save(strcat('delta_waveform_broadband_',g{j},'.mat'),'delta_waveform_broadband_total','-v7.3')
  save(strcat('delta_waveform_',g{j},'.mat'),'delta_waveform_total','-v7.3')
  delta_timestamps_SD =  delta_total_data{j}.';
  d_count = cellfun(@size,delta_timestamps_SD,'UniformOutput',false);
  for d = 1:length(d_count)
      delta_count(d) = d_count{d}(1);
      if isnan(delta_timestamps_SD{d})
          delta_count(d) = 0;
      end 
  end 
  save(strcat('delta_timestamps_',g{j},'.mat'),'delta_timestamps_SD','-v7.3')
  save(strcat('delta_count_',g{j},'.mat'),'delta_count','-v7.3')
  delta_waveform_broadband_comp = [delta_waveform_broadband_comp;delta_waveform_broadband_total];
  delta_waveform_comp = [delta_waveform_comp;delta_waveform_total];
  save(strcat('delta_waveform_broadband_compilation_Rat',rat_folder{k},'.mat'),'delta_waveform_broadband_comp','-v7.3')
  save(strcat('delta_waveform_compilation_Rat',rat_folder{k},'.mat'),'delta_waveform_comp','-v7.3')
    end
    cd ..
end

        
%     if ind_mode==1
%         %% Threshold selection
%         [swr_hpc,swr_pfc,s_hpc,s_pfc,V_hpc,V_pfc,signal2_hpc,signal2_pfc,sd_swr]=swr_check_thr(HPC,Cortex,states,ss,D1,D2);       
%         max_length=cellfun(@length,V_hpc);
%         N=max_length==max(max_length);
%         % max_length=cellfun(@length,swr_pfc(:,1));
%         % N=max_length==max(max_length);
% 
%         hpc=V_hpc{N};
%         pfc=V_pfc{N};
%         hpc2=signal2_hpc{N};
%         pfc2=signal2_pfc{N};
%         n=find(N);
% 
%         plot((1:length(hpc))./1000,0.1.*zscore(hpc)+100,'Color','black')
%         hold on
%         plot((1:length(pfc))./1000,0.1.*zscore(pfc)+150,'Color','black')
%         xlabel('Time (Seconds)')
% 
% 
%         stem([swr_hpc{n,3}],ones(length([swr_hpc{n}]),1).*350,'Color','blue') %(HPC)
%         stem([swr_pfc{n,3}],ones(length([swr_pfc{n}]),1).*350,'Color','red')%Seconds (Cortex)
% 
%         % 
%         plot((1:length(hpc2))./1000,0.1.*zscore(hpc2)+220,'Color','black')
%         plot((1:length(pfc2))./1000,0.1.*zscore(pfc2)+290,'Color','black')
%         % 
%         yticks([100 150 220 290])
%         yticklabels({'HPC',xx{1},'HPC (Bandpassed)',[xx{1}
%         '(Bandpassed)']})k
%         a = get(gca,'YTickLabel');
%         set(gca,'YTickLabel',a,'FontName','Times','fontsize',12)
%     'Stop the code here'
%     
%             %% Threshold selection
%         [nr_swr_HPC, nr_swr_Cortex,nr_cohfos,nr_single_hpc,nr_single_cortex,place1,place2,V_hpc,V_pfc,signal2_hpc,signal2_pfc,HPC_binned,Cortex_binned]=bin_swr_detection(HPC,Cortex,states,ss,D1,D2,xx,yy,fn);       
% 
%         % max_length=cellfun(@length,swr_pfc(:,1));
%         % N=max_length==max(max_length);
% 
%         hpc=V_hpc{place1}{place2};
%         pfc=V_pfc{place1}{place2};
%         hpc2=signal2_hpc{place1}{place2};
%         pfc2=signal2_pfc{place1}{place2};
% %          n=find(N);
% 
%         plot((1:length(hpc))./1000,0.1.*zscore(hpc)+100,'Color','black')
%         hold on
%         plot((1:length(pfc))./1000,0.1.*zscore(pfc)+150,'Color','black')
%         xlabel('Time (Seconds)')
% 
% 
%         stem([HPC_binned{place1}{place2,3}].',ones(length([HPC_binned{place1}{place2}]),1).*350,'Color','blue') %(HPC)
%         stem([Cortex_binned{place1}{place2,3}].',ones(length([Cortex_binned{place1}{place2}]),1).*350,'Color','red')%Seconds (Cortex)
% 
%         % 
%         plot((1:length(hpc2))./1000,0.1.*zscore(hpc2)+220,'Color','black')
%         plot((1:length(pfc2))./1000,0.1.*zscore(pfc2)+290,'Color','black')
%         % 
%         yticks([100 150 220 290])
%         yticklabels({'HPC',xx{1},'HPC (Bandpassed)',[xx{1} '(Bandpassed)']})
%         a = get(gca,'YTickLabel');
%         set(gca,'YTickLabel',a,'FontName','Times','fontsize',12)
%     'Stop the code here'
%         xo
%     end
    
    
%     if ind_mode==2
%         %Count events
%         [nr_swr_HPC(i,:), nr_swr_Cortex(i,:),nr_cohfos(i,:)]=bin_swr_detection(HPC,Cortex,states,ss,D1,D2,xx,yy,fn); 
%     end

%     
%     if ind_mode==3
%         %Find SD values
%         [sd_swr]=find_std(HPC,Cortex,states,ss);
%         
%         Sd_Swr.sd2_hpc_co(i)=sd_swr.sd2_hpc_co;
%         Sd_Swr.sd5_hpc_co(i)=sd_swr.sd5_hpc_co;
%         Sd_Swr.sd2_pfc_co(i)=sd_swr.sd2_pfc_co;
%         Sd_Swr.sd5_pfc_co(i)=sd_swr.sd5_pfc_co;
%         Sd_Swr.sd2_hpc_long(i)=sd_swr.sd2_hpc_long;
%         Sd_Swr.sd5_hpc_long(i)=sd_swr.sd5_hpc_long;
%         Sd_Swr.sd2_pfc_long(i)=sd_swr.sd2_pfc_long;
%         Sd_Swr.sd5_pfc_long(i)=sd_swr.sd5_pfc_long;
%         
% %         if i==5
% %             xo
% %         end
%     end
%               
%     
%     
%                     elseif contains(G{i}, 'rial5') % PostTrial 5 case 
% %                         xo
%                         %Sleep scoring data
%                         if length(states)<45*60*4
%                             states=[states nan(1,45*60*4-length(states))]; %Fill with NaNs.
%                         else
%                             states=states(1:45*60*4); %Take only 45 min.
%                         end
%                         
%                         
%                         %Ephys
%                         if length(HPC)<45*60*1000*4
%                             HPC=[HPC.' (nan(45*60*1000*4-length(HPC),1).')]; %Fill with NaNs.
%                         else
%                             HPC=HPC(1:45*60*1000*4).'; %Take only 45 min.
%                         end
%                         
%                         if length(Cortex)<45*60*1000*4
%                             Cortex=[Cortex.' (nan(45*60*1000*4-length(Cortex),1).')]; %Fill with NaNs.
%                         else
%                             Cortex=Cortex(1:45*60*1000*4).'; %Take only 45 min.
%                         end
% 
%                         
%                         %Chunk in 4
%                         states1=states(1:2700);
%                         states2=states(2700+1:2700*2);
%                         states3=states(1+2700*2:2700*3);
%                         states4=states(1+2700*3:2700*4);
%                         
%                         HPC_1=HPC(1:2700*1000);
%                         HPC_2=HPC(2700*1000+1:2700*2*1000);
%                         HPC_3=HPC(1+2700*2*1000:2700*3*1000);
%                         HPC_4=HPC(1+2700*3*1000:2700*4*1000);
%                         
%                         Cortex_1=Cortex(1:2700*1000);
%                         Cortex_2=Cortex(2700*1000+1:2700*2*1000);
%                         Cortex_3=Cortex(1+2700*2*1000:2700*3*1000);
%                         Cortex_4=Cortex(1+2700*3*1000:2700*4*1000);
%                         
%                         if ind_mode==2
%                             [nr_swr_HPC(6,:), nr_swr_Cortex(6,:),nr_cohfos(6,:)]=bin_swr_detection(HPC_1,Cortex_1,states1,ss,D1,D2,xx,yy,fn);  
%                             [nr_swr_HPC(7,:), nr_swr_Cortex(7,:),nr_cohfos(7,:)]=bin_swr_detection(HPC_2,Cortex_2,states2,ss,D1,D2,xx,yy,fn);  
%                             [nr_swr_HPC(8,:), nr_swr_Cortex(8,:),nr_cohfos(8,:)]=bin_swr_detection(HPC_3,Cortex_3,states3,ss,D1,D2,xx,yy,fn);  
%                             [nr_swr_HPC(9,:), nr_swr_Cortex(9,:),nr_cohfos(9,:)]=bin_swr_detection(HPC_4,Cortex_4,states4,ss,D1,D2,xx,yy,fn);  
%                         end
%                         
%                         
%                         if ind_mode==3
%                                     %Find SD values
%                                     %[~,~,~,~,~,~,~,~,sd_swr]=swr_check_thr(HPC_1,Cortex_1,states1,ss,70,30);
%                                     [sd_swr]=find_std(HPC_1,Cortex_1,states1,ss);
%                                     Sd_Swr.sd2_hpc_co(6)=sd_swr.sd2_hpc_co;
%                                     Sd_Swr.sd5_hpc_co(6)=sd_swr.sd5_hpc_co;
%                                     Sd_Swr.sd2_pfc_co(6)=sd_swr.sd2_pfc_co;
%                                     Sd_Swr.sd5_pfc_co(6)=sd_swr.sd5_pfc_co;
%                                     Sd_Swr.sd2_hpc_long(6)=sd_swr.sd2_hpc_long;
%                                     Sd_Swr.sd5_hpc_long(6)=sd_swr.sd5_hpc_long;
%                                     Sd_Swr.sd2_pfc_long(6)=sd_swr.sd2_pfc_long;
%                                     Sd_Swr.sd5_pfc_long(6)=sd_swr.sd5_pfc_long;
%                                     
%                                     
%                                     %[~,~,~,~,~,~,~,~,sd_swr]=swr_check_thr(HPC_2,Cortex_2,states2,ss,70,30);
%                                     [sd_swr]=find_std(HPC_2,Cortex_2,states2,ss);
%                                     Sd_Swr.sd2_hpc_co(7)=sd_swr.sd2_hpc_co;
%                                     Sd_Swr.sd5_hpc_co(7)=sd_swr.sd5_hpc_co;
%                                     Sd_Swr.sd2_pfc_co(7)=sd_swr.sd2_pfc_co;
%                                     Sd_Swr.sd5_pfc_co(7)=sd_swr.sd5_pfc_co;
%                                     Sd_Swr.sd2_hpc_long(7)=sd_swr.sd2_hpc_long;
%                                     Sd_Swr.sd5_hpc_long(7)=sd_swr.sd5_hpc_long;
%                                     Sd_Swr.sd2_pfc_long(7)=sd_swr.sd2_pfc_long;
%                                     Sd_Swr.sd5_pfc_long(7)=sd_swr.sd5_pfc_long;
%                                     
%                                     %[~,~,~,~,~,~,~,~,sd_swr]=swr_check_thr(HPC_3,Cortex_3,states3,ss,70,30);
%                                     [sd_swr]=find_std(HPC_3,Cortex_3,states3,ss);
%                                     Sd_Swr.sd2_hpc_co(8)=sd_swr.sd2_hpc_co;
%                                     Sd_Swr.sd5_hpc_co(8)=sd_swr.sd5_hpc_co;
%                                     Sd_Swr.sd2_pfc_co(8)=sd_swr.sd2_pfc_co;
%                                     Sd_Swr.sd5_pfc_co(8)=sd_swr.sd5_pfc_co;
%                                     Sd_Swr.sd2_hpc_long(8)=sd_swr.sd2_hpc_long;
%                                     Sd_Swr.sd5_hpc_long(8)=sd_swr.sd5_hpc_long;
%                                     Sd_Swr.sd2_pfc_long(8)=sd_swr.sd2_pfc_long;
%                                     Sd_Swr.sd5_pfc_long(8)=sd_swr.sd5_pfc_long;       
%                                     
%                                     %[~,~,~,~,~,~,~,~,sd_swr]=swr_check_thr(HPC_4,Cortex_4,states4,ss,70,30);
%                                     [sd_swr]=find_std(HPC_4,Cortex_4,states4,ss);
%                                     Sd_Swr.sd2_hpc_co(9)=sd_swr.sd2_hpc_co;
%                                     Sd_Swr.sd5_hpc_co(9)=sd_swr.sd5_hpc_co;
%                                     Sd_Swr.sd2_pfc_co(9)=sd_swr.sd2_pfc_co;
%                                     Sd_Swr.sd5_pfc_co(9)=sd_swr.sd5_pfc_co;
%                                     Sd_Swr.sd2_hpc_long(9)=sd_swr.sd2_hpc_long;
%                                     Sd_Swr.sd5_hpc_long(9)=sd_swr.sd5_hpc_long;
%                                     Sd_Swr.sd2_pfc_long(9)=sd_swr.sd2_pfc_long;
%                                     Sd_Swr.sd5_pfc_long(9)=sd_swr.sd5_pfc_long;                                      
% 
%                         end
% 
%     
%                     end
% 
% 
% 
% %%
% 
%                     cd ..
%                 else
%                     cd .. %Means there is no sleep scoring file.
%                     end
%                 else
%                     cd ..
%             end
% % if i==5
% %     xo
% % end            
%             end
%                 cd ..
% 
%         end
%         
%    if ind_mode==2
%     for d=1:size(nr_swr_Cortex,1)
%         SD_nr_Cortex(d)=sum(nr_swr_Cortex(d,:));
%         SD_nr_HPC(d)=sum(nr_swr_HPC(d,:));
%         SD_nr_cohfos(d)=sum(nr_cohfos(d,:));
%     end
% 
% xo
% 
% %         allscreen()
% %         subplot(1,2,1)
%         imagesc(nr_swr_HPC);  olorbar(); colormap('gray')
%         xticks([1.5:9.5])
%         xticklabels({'5','10','15','20','25','30','35','40','45'})
%         yticks([1:9])
%         yticklabels({'PS','PT1','PT2','PT3','PT4','PT5_1','PT5_2','PT5_3','PT5_4'})
%         title('HPC HFOs')
%         printing([yy{1} '_' g{j} '_' stage])
%         close all
%         
% %         subplot(1,2,2)
%         imagesc(nr_swr_Cortex); colorbar(); colormap('gray')
%         xticks([1.5:9.5])
%         xticklabels({'5','10','15','20','25','30','35','40','45'})
%         yticks([1:9])
%         yticklabels({'PS','PT1','PT2','PT3','PT4','PT5_1','PT5_2','PT5_3','PT5_4'})
%         title([xx{1} ' HFOs'])
%         printing([xx{1} '_' g{j} '_' stage])
%         close all       
%    end
%     
%    if ind_mode==3
% %        xo
%        %Sd_Swr
%     TT=table;
%     TT.Variables=...
%         [
%                 [{g{j}};{'x'};{'x'};{'x'};{'x'};{'x'};{'x'};{'x'}] ...
%         [{'HPC_2SD_Concatenated'};{'HPC_2SD_Longest'};{'HPC_5SD_Concatenated'};{'HPC_5SD_Longest'};{'PFC_2SD_Concatenated'};{'PFC_2SD_Longest'};{'PFC_5SD_Concatenated'};{'PFC_5SD_Longest'}] ...
%     num2cell([ Sd_Swr.sd2_hpc_co;Sd_Swr.sd2_hpc_long;Sd_Swr.sd5_hpc_co;Sd_Swr.sd5_hpc_long ...
%         ;Sd_Swr.sd2_pfc_co;Sd_Swr.sd2_pfc_long;Sd_Swr.sd5_pfc_co;Sd_Swr.sd5_pfc_long ...
%       ]) ...
%     num2cell([ mean(Sd_Swr.sd2_hpc_co); mean(Sd_Swr.sd2_hpc_long); mean(Sd_Swr.sd5_hpc_co);mean(Sd_Swr.sd5_hpc_long) ...
%         ;mean(Sd_Swr.sd2_pfc_co);mean(Sd_Swr.sd2_pfc_long);mean(Sd_Swr.sd5_pfc_co);mean(Sd_Swr.sd5_pfc_long) ...
%       ])...
%       ];
% corrected_means = []; 
%  TT2=cell2mat(TT{:,3:end-1});
%  for m=1:size(TT2,1)
%     corrected_means = [corrected_means; mean(rmoutliers(nonzeros(TT2(m,3:end))))];
%  end
%  TT.corrected_means = num2cell(corrected_means);
%  
% %     TT.Properties.VariableNames={'Condition';'Metric';'PS';'PT1';'PT2';'PT3';'PT4';'PT5_1';'PT5_2';'PT5_3';'PT5_4';'Average_Value';'Average_Outliers_Removed'};
% %     t1=repmat({'x'},[1 12]);
%     
% %         if j==1
% %             T=[TT];
% %         else
% %             T=[];
% %             T=[T;t1;TT];
% %         end
% %     
% 
%        cd(g{j})
%             for i=1:length(G);
% %             for i=1
% %                  xo
%                 cd(G{i})
%                 clear states
%                 clear HPC Cortex
%                 A = dir('*states*.mat');
%                 A={A.name};
%                 
%                 if sum(contains(A, 'states')) > 0 %More than 2 sleep scoring files
%                     A=A(cellfun(@(x) ~isempty(strfind(x,'states')),A));
%                     A=A(~(cellfun(@(x) ~isempty(strfind(x,'eeg')),A)));
%                     if sum(contains(A, 'states')) > 0
%                     
% %                     st2=st(cellfun(@(x) ~isempty(strfind(x,barea)),st)); %Brain area.
%                     cellfun(@load,A);
% 
% 
%                    
%                     HPC=dir('*HPC_*.mat');
%                     HPC=HPC.name;
%                     HPC=load(HPC);
%                     HPC=HPC.HPC;
%                     HPC=HPC.*(0.195);
%                     
%                     Cortex=dir(strcat('*',xx{1},'*.mat'));
%                     Cortex=Cortex.name;
%                     Cortex=load(Cortex);
%                     % Cortex=Cortex.Cortex;
%                     Cortex=getfield(Cortex,xx{1});
%                     Cortex=Cortex.*(0.195);
% 
% 
% 
%                                       
% 
%                     if and(~contains(G{i},'trial5'),~contains(G{i},'Trial5')) %Whenever it is not PostTrial 5 
%                         
%                         % Sleep Scoring data
%                         if length(states)<45*60
%                             states=[states nan(1,45*60-length(states))]; %Fill with NaNs.
%                         else
%                             states=states(1:45*60); %Take only 45 min.
%                         end
%                         
%                         %Ephys data
%                         if length(HPC)<45*60*1000
%                             HPC=[HPC.' (nan(45*60*1000-length(HPC),1).')]; %Fill with NaNs.
%                         else
%                             HPC=HPC(1:45*60*1000).'; %Take only 45 min.
%                         end
%                         
%                         if length(Cortex)<45*60*1000
%                             Cortex=[Cortex.' (nan(45*60*1000-length(Cortex),1).')]; %Fill with NaNs.
%                         else
%                             Cortex=Cortex(1:45*60*1000).'; %Take only 45 min.
%                         end
% 
%                         D1 = round(cell2mat(TT{3,end}) + str2num(cell2mat(offset1)));
%                         D2 = round(cell2mat(TT{7,end}) + str2num(cell2mat(offset2)));
%                         [swr_hpc,swr_pfc,s_hpc,s_pfc,V_hpc,V_pfc,signal2_hpc,signal2_pfc,sd_swr, M_multiplets, Mx_multiplets, multiplets]=swr_check_thr(HPC,Cortex,states,ss,D1,D2);
% %                         [nr_swr_HPC(i,:), nr_swr_Cortex(i,:),nr_cohfos(i,:)]=bin_swr_detection(HPC,Cortex,states,ss,D1,D2,xx,yy,fn);
%                             multiplet_count={};
%                             for l=1:length(multiplets)
%                                 multiplet_count{l}=sum(cellfun('length', Mx_multiplets.(multiplets{l}).m_1));
% 
%                             end
%                             total_swrs(j,i) = sum(s_hpc);
%                             total_hfos(j,i) = sum(s_pfc);
%                             multiplet_count_total{j,i}=multiplet_count;
%                             stage_count = sum(states(:)==ss);
%                             total_swrs_minute(j,i)=(total_swrs(j,i)/stage_count*60);
%                             total_hfos_minute(j,i)=(total_hfos(j,i)/stage_count*60);
%                             
% %                                     [nr_swr_HPC(i,:), nr_swr_Cortex(i,:), nr_cohfos(i,:),nr_single_hpc(i,:),nr_single_cortex(i,:),mfreq_hpc,mfreq_cortex]=bin_swr_detection(HPC,Cortex,states,ss,D1,D2,xx,yy,fn); 
% %         mfreq_hpc=mfreq_hpc(~cellfun('isempty',mfreq_hpc));
% %         mfreq_cortex=mfreq_cortex(~cellfun('isempty',mfreq_cortex));
% %         mfreq_hpc=[mfreq_hpc{:}];
% %         mfreq_cortex=[mfreq_cortex{:}];
% %         Mfreq_hpc{i}=mfreq_hpc;
% %         Mfreq_cortex{i}=mfreq_cortex;
% %                             
%                             
% 
%                             end
% 
%                                                 if contains(G{i}, 'rial5')
%                                                     clear HPC_1 HPC_2 HPC_3 HPC_4 Cortex_1 Cortex_2 Cortex_3 Cortex_4
%                                                     
%                                                                             %Sleep scoring data
%                         if length(states)<45*60*4
%                             states=[states nan(1,45*60*4-length(states))]; %Fill with NaNs.
%                         else
%                             states=states(1:45*60*4); %Take only 45 min.
%                         end
%                         
%                         
%                         %Ephys
%                         if length(HPC)<45*60*1000*4
%                             HPC=[HPC.' (nan(45*60*1000*4-length(HPC),1).')]; %Fill with NaNs.
%                         else
%                             HPC=HPC(1:45*60*1000*4).'; %Take only 45 min.
%                         end
%                         
%                         if length(Cortex)<45*60*1000*4
%                             Cortex=[Cortex.' (nan(45*60*1000*4-length(Cortex),1).')]; %Fill with NaNs.
%                         else
%                             Cortex=Cortex(1:45*60*1000*4).'; %Take only 45 min.
%                         end
% 
%                         
%                         %Chunk in 4
%                         states1=states(1:2700);
%                         states2=states(2700+1:2700*2);
%                         states3=states(1+2700*2:2700*3);
%                         states4=states(1+2700*3:2700*4);
%                         
%                         HPC_1=HPC(1:2700*1000);
%                         HPC_2=HPC(2700*1000+1:2700*2*1000);
%                         HPC_3=HPC(1+2700*2*1000:2700*3*1000);
%                         HPC_4=HPC(1+2700*3*1000:2700*4*1000);
%                         
%                         Cortex_1=Cortex(1:2700*1000);
%                         Cortex_2=Cortex(2700*1000+1:2700*2*1000);
%                         Cortex_3=Cortex(1+2700*2*1000:2700*3*1000);
%                         Cortex_4=Cortex(1+2700*3*1000:2700*4*1000);
%                             
%                                 [swr_hpc,swr_pfc,s_hpc,s_pfc,V_hpc,V_pfc,signal2_hpc,signal2_pfc,sd_swr, M_multiplets, Mx_multiplets, multiplets]=swr_check_thr(HPC_1,Cortex_1,states1,ss,D1,D2);
% %                         [nr_swr_HPC(i,:), nr_swr_Cortex(i,:),nr_cohfos(i,:)]=bin_swr_detection(HPC,Cortex,states,ss,D1,D2,xx,yy,fn);
%                                 multiplet_count={};
%                                 for l=1:length(multiplets)
%                                     multiplet_count{l}=sum(cellfun('length', Mx_multiplets.(multiplets{l}).m_1));
% 
%                                 end
%                                 
%                                 total_swrs(j,i) = sum(s_hpc);
%                                 total_hfos(j,i) = sum(s_pfc);
%                                 stage_count = sum(states1(:)==ss);
%                                 multiplet_count_total{j,i}=multiplet_count;
%                                 total_swrs_minute(j,i)=(total_swrs(j,i)/stage_count*60);
%                                 total_hfos_minute(j,i)=(total_hfos(j,i)/stage_count*60);
% 
%                                 [swr_hpc,swr_pfc,s_hpc,s_pfc,V_hpc,V_pfc,signal2_hpc,signal2_pfc,sd_swr, M_multiplets, Mx_multiplets, multiplets]=swr_check_thr(HPC_2,Cortex_2,states2,ss,D1,D2);
%                                 multiplet_count={};
%                                 for l=1:length(multiplets)
%                                     multiplet_count{l}=sum(cellfun('length', Mx_multiplets.(multiplets{l}).m_1));
% 
%                                 end
%                                 
%                                 total_swrs(j,i+1) = sum(s_hpc);
%                                 total_hfos(j,i+1) = sum(s_pfc);
%                                 stage_count = sum(states2(:)==ss);
%                                 multiplet_count_total{j,i+1}=multiplet_count;
%                                 total_swrs_minute(j,i+1)=(total_swrs(j,i+1)/stage_count*60);
%                                 total_hfos_minute(j,i+1)=(total_hfos(j,i+1)/stage_count*60);
% 
%                                 [swr_hpc,swr_pfc,s_hpc,s_pfc,V_hpc,V_pfc,signal2_hpc,signal2_pfc,sd_swr, M_multiplets, Mx_multiplets, multiplets]=swr_check_thr(HPC_3,Cortex_3,states3,ss,D1,D2);
%                                 multiplet_count={};
%                                 for l=1:length(multiplets)
%                                     multiplet_count{l}=sum(cellfun('length', Mx_multiplets.(multiplets{l}).m_1));
% 
%                                 end
%                                 
%                                 total_swrs(j,i+2) = sum(s_hpc);
%                                 total_hfos(j,i+2) = sum(s_pfc);
%                                 stage_count = sum(states3(:)==ss);
%                                 multiplet_count_total{j,i+2}=multiplet_count;
%                                 total_swrs_minute(j,i+2)=(total_swrs(j,i+2)/stage_count*60);
%                                 total_hfos_minute(j,i+2)=(total_hfos(j,i+2)/stage_count*60);
% 
%                                 [swr_hpc,swr_pfc,s_hpc,s_pfc,V_hpc,V_pfc,signal2_hpc,signal2_pfc,sd_swr, M_multiplets, Mx_multiplets, multiplets]=swr_check_thr(HPC_4,Cortex_4,states4,ss,D1,D2);
%                                 multiplet_count={};
%                                 for l=1:length(multiplets)
%                                     multiplet_count{l}=sum(cellfun('length', Mx_multiplets.(multiplets{l}).m_1));
% 
%                                 end
%                                 
%                                 total_swrs(j,i+3) = sum(s_hpc);
%                                 total_hfos(j,i+3) = sum(s_pfc);
%                                 stage_count = sum(states4(:)==ss);
%                                 multiplet_count_total{j,i+3}=multiplet_count;
%                                 total_swrs_minute(j,i+3)=(total_swrs(j,i+3)/stage_count*60);
%                                 total_hfos_minute(j,i+3)=(total_hfos(j,i+3)/stage_count*60);
%                                 
% %                             [nr_swr_HPC(6,:), nr_swr_Cortex(6,:),nr_cohfos(6,:)]=bin_swr_detection(HPC_1,Cortex_1,states1,ss,D1,D2,xx,yy,fn);  
% %                             [nr_swr_HPC(7,:), nr_swr_Cortex(7,:),nr_cohfos(7,:)]=bin_swr_detection(HPC_2,Cortex_2,states2,ss,D1,D2,xx,yy,fn);  
% %                             [nr_swr_HPC(8,:), nr_swr_Cortex(8,:),nr_cohfos(8,:)]=bin_swr_detection(HPC_3,Cortex_3,states3,ss,D1,D2,xx,yy,fn);  
% %                             [nr_swr_HPC(9,:), nr_swr_Cortex(9,:),nr_cohfos(9,:)]=bin_swr_detection(HPC_4,Cortex_4,states4,ss,D1,D2,xx,yy,fn);
% % 
% %                                                         [nr_swr_HPC(6,:), nr_swr_Cortex(6,:),nr_cohfos(6,:),nr_single_hpc(6,:),nr_single_cortex(6,:),mfreq_hpc,mfreq_cortex]=bin_swr_detection(HPC_1,Cortex_1,states1,ss,D1,D2,xx,yy,fn);
% %                                 mfreq_hpc=mfreq_hpc(~cellfun('isempty',mfreq_hpc));
% %                                 mfreq_cortex=mfreq_cortex(~cellfun('isempty',mfreq_cortex));
% %                                 mfreq_hpc=[mfreq_hpc{:}];
% %                                 mfreq_cortex=[mfreq_cortex{:}];
% %                                 Mfreq_hpc{6}=mfreq_hpc;
% %                                 Mfreq_cortex{6}=mfreq_cortex;  
% % 
% %                             [nr_swr_HPC(7,:), nr_swr_Cortex(7,:),nr_cohfos(7,:),nr_single_hpc(7,:),nr_single_cortex(7,:),mfreq_hpc,mfreq_cortex]=bin_swr_detection(HPC_2,Cortex_2,states2,ss,D1,D2,xx,yy,fn);
% %                                 mfreq_hpc=mfreq_hpc(~cellfun('isempty',mfreq_hpc));
% %                                 mfreq_cortex=mfreq_cortex(~cellfun('isempty',mfreq_cortex));
% %                                 mfreq_hpc=[mfreq_hpc{:}];
% %                                 mfreq_cortex=[mfreq_cortex{:}];
% %                                 Mfreq_hpc{7}=mfreq_hpc;
% %                                 Mfreq_cortex{7}=mfreq_cortex;  
% % 
% %                             [nr_swr_HPC(8,:), nr_swr_Cortex(8,:),nr_cohfos(8,:),nr_single_hpc(8,:),nr_single_cortex(8,:),mfreq_hpc,mfreq_cortex]=bin_swr_detection(HPC_3,Cortex_3,states3,ss,D1,D2,xx,yy,fn);  
% %                                 mfreq_hpc=mfreq_hpc(~cellfun('isempty',mfreq_hpc));
% %                                 mfreq_cortex=mfreq_cortex(~cellfun('isempty',mfreq_cortex));
% %                                 mfreq_hpc=[mfreq_hpc{:}];
% %                                 mfreq_cortex=[mfreq_cortex{:}];
% %                                 Mfreq_hpc{8}=mfreq_hpc;
% %                                 Mfreq_cortex{8}=mfreq_cortex;
% % 
% %                             [nr_swr_HPC(9,:), nr_swr_Cortex(9,:),nr_cohfos(9,:),nr_single_hpc(9,:),nr_single_cortex(9,:),mfreq_hpc,mfreq_cortex]=bin_swr_detection(HPC_4,Cortex_4,states4,ss,D1,D2,xx,yy,fn);
% %                                 mfreq_hpc=mfreq_hpc(~cellfun('isempty',mfreq_hpc));
% %                                 mfreq_cortex=mfreq_cortex(~cellfun('isempty',mfreq_cortex));
% %                                 mfreq_hpc=[mfreq_hpc{:}];
% %                                 mfreq_cortex=[mfreq_cortex{:}];
% %                                 Mfreq_hpc{9}=mfreq_hpc;
% %                                 Mfreq_cortex{9}=mfreq_cortex;                      
% %                             
%                                                 end
% 
%                     end
%                 end
%                     
% 
%                         
%                 cd ..
%                 end
%                 cd ..
%                 D_all(1,j) = D1;
%                 D_all(2,j) = D2;
% %                 nr_cohfos_pt_animal(j,:)=nr_cohfos_pt(1,:);
%                        if ind_mode==3
% % %xo
% % if k==3 || k==7
% %     sdpos=strfind(g{j},'SD');
% %     if (contains(g{j}, 'OR') && ~contains(g{j}, 'NOV'))
% %     condition = 'OR';
% %     end
% %     if contains(g{j}, 'OD')
% %         condition = 'OD';
% %     end
% %     if contains(g{j}, 'CON')
% %         condition = 'CON';
% %     end
% %     if contains(g{j}, 'HC')
% %         condition = 'HC';
% %     end
% %     if contains(g{j}, 'OR_N')
% %         condition = 'OR_N';
% %     end
% % 
% %     save(['Mfreq_' xx{1} '_' condition '_' g{j}(sdpos:sdpos+3) '.mat'],'Mfreq_cortex');
% %     save(['Mfreq_' yy{1} '_' condition '_' g{j}(sdpos:sdpos+3) '.mat'],'Mfreq_hpc');
% % 
% % else
% %     save(['Mfreq_' xx{1} '_' g{j} '.mat'],'Mfreq_cortex');
% %     save(['Mfreq_' yy{1} '_' g{j} '.mat'],'Mfreq_hpc');
% % end
% % 
% % clear Mfreq_cortex Mfreq_hpc
%     
% 
% 
%     end
%    end
% %            SD_nr_Cortex=zeros(1,9);
% %         SD_nr_HPC=zeros(1,9);
% %         SD_nr_cohfos=zeros(1,9);
% %         for d=1:size(nr_swr_Cortex,1)
% %             SD_nr_Cortex(d)=sum(nr_swr_Cortex(d,:));
% %             SD_nr_HPC(d)=sum(nr_swr_HPC(d,:));
% %             SD_nr_cohfos(d)=sum(nr_cohfos(d,:));
% %         end
% %             animal_nr_Cortex(j,:)=SD_nr_Cortex;
% %             animal_nr_HPC(j,:)=SD_nr_HPC;
% %             animal_nr_cohfos(j,:)=SD_nr_cohfos;
%             
%             
%               
% %           clear SD_nr_Cortex SD_nr_HPC SD_nr_cohfos nr_swr_HPC nr_swr_Cortex nr_cohfos
%     end
% %        if ind_mode==3 
% %            if k==1
% %                 P_bar=waitbar(0,'Please wait...');
% %            end
% %                progress_bar(k,length(rat_folder),P_bar)
% %         end  
%     %xo
%      %writetable(T,strcat('SD_values','.xls'),'Sheet',1,'Range','A1:Z50')
%     cd ..
% end
% 
%         zsinglet_total=[];
%         zdoublet_total=[];
%         ztriplet_total=[];
%         zquadruplet_total=[];
%         zquintuplet_total=[];
%         zsextuplet_total=[];
%         zseptuplet_total=[];
%         zoctuplet_total=[];
%         znonuplet_total=[];
% 
% for r=1:size(multiplet_count_total,2)
%     for e=1:size(multiplet_count_total,1)
%         if ~isempty(multiplet_count_total{e,r})
%             zsinglet_total(e,r)=cell2mat(multiplet_count_total{e,r}(1,1));
%             zdoublet_total(e,r)=cell2mat(multiplet_count_total{e,r}(1,2));
%             ztriplet_total(e,r)=cell2mat(multiplet_count_total{e,r}(1,3));
%             zquadruplet_total(e,r)=cell2mat(multiplet_count_total{e,r}(1,4));
%             zquintuplet_total(e,r)=cell2mat(multiplet_count_total{e,r}(1,5));
%             zsextuplet_total(e,r)=cell2mat(multiplet_count_total{e,r}(1,6));
%             zseptuplet_total(e,r)=cell2mat(multiplet_count_total{e,r}(1,7));
%             zoctuplet_total(e,r)=cell2mat(multiplet_count_total{e,r}(1,8));
%             znonuplet_total(e,r)=cell2mat(multiplet_count_total{e,r}(1,9));
%         end
%     end
% end



