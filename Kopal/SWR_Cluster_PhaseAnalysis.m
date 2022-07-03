clear variables
addpath(genpath('/Volumes/Samsung_T5/Milan_DA/OS_ephys_da/CorticoHippocampal'))
addpath ('/Volumes/Samsung_T5/Milan_DA/OS_ephys_da/ADRITOOLS')
 cd('/Volumes/Samsung_T5/Milan_DA/RGS14_Ephys_da/Data_RGS14_Downsampled_First_Session')
      yy = {'HPC'};       
      xx = {'PFC'};
      ss = 3;   %NREM
%% Threshold Offset
prompt = {'Select a threshold offset from 5SD for HPC'};
dlgtitle = 'Threshold offset HPC';
definput = {'5'};
offset1 = inputdlg(prompt,dlgtitle,[1 40],definput);

prompt = {'Select a threshold offset from 5SD for Cortex'};
dlgtitle = 'Threshold Cortex';
definput = {'5'};
offset2 = inputdlg(prompt,dlgtitle,[1 40],definput);

nr_swr_HPC = [];
nr_swr_Cortex = [];
D_all = [];
nr_cohfos_pt_animal = [];    
fn = 1000;
rat_folder = getfolder;
nr_cohfos_pt = [];
% for k=1:length(rat_folder)
for k = 4
     %rat index 
  cd(rat_folder{k})
ripple_phase_comp_c1 = [];
ripple_phase_comp_c2 = [];
ripple_phase_comp_c3 = [];
ripple_phase_comp_c4 = [];

g = getfolder;
    for j = 1:length(g)
    %% Finding SD data     
        current_sd = g{j};
        current_sd = strsplit(current_sd, '_');
        idx = find(contains(current_sd, 'sd', 'IgnoreCase',true));
        SDn = current_sd{idx}(1,3:end);
        cd ('/Volumes/Samsung_T5/Milan_DA/SWR_Cluster_SDwise_pelin')
        ClSDw =  getfolder;
        cd(ClSDw{k})
        path = cd;
        dinfo = dir(path);
        dinfo = {dinfo.name};
        dinfo = dinfo(3:end);
        idx2 = find(contains(dinfo, ['SD',SDn,'_'], 'IgnoreCase',true));
        dinfo_sd = dinfo(idx2);
       clus1 = load(dinfo_sd{contains(dinfo_sd,'c1')});
       name = fieldnames(clus1);
       clus1 = clus1.(name{1});
       clus2 = load(dinfo_sd{contains(dinfo_sd,'c2')});
       name = fieldnames(clus2);
       clus2 = clus2.(name{1});
       clus3 = load(dinfo_sd{contains(dinfo_sd,'c3')});
        name = fieldnames(clus3);
       clus3 = clus3.(name{1});
       clus4 = load(dinfo_sd{contains(dinfo_sd,'c4')});
        name = fieldnames(clus4);
       clus4 = clus4.(name{1});
            clus1(:,1) = [];
            clus2(:,1) = [];
            clus3(:,1) = [];
            clus4(:,1) = [];
 %% Loading Ephys data
 cd('/Volumes/Samsung_T5/Milan_DA/RGS14_Ephys_da/Data_RGS14_Downsampled_First_Session') 
 rat_folder = getfolder;
 cd(rat_folder{k})  
        
 cd(g{j})
 G=getfolder;
 r_cohfos_pt=zeros(1,9);              
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

    for i=1:length(G); 
%              for i = 3 
        clear states
%
%                 for i=3
%                  xo
        cd(G{i})
        A = dir('*states*.mat');
        A={A.name};

        if sum(contains(A, 'states')) > 0 %More than 2 sleep scoring files
            A=A(cellfun(@(x) ~isempty(strfind(x,'states')),A));
            A=A(~(cellfun(@(x) ~isempty(strfind(x,'eeg')),A)));

            if sum(contains(A, 'states')) > 0
            cellfun(@load,A);


HPC=dir('*HPC_*.mat');
HPC=HPC.name;
HPC=load(HPC);
HPC=HPC.HPC;
HPC=HPC.*(0.195);

Cortex=dir(strcat('*',xx{1},'*.mat'));
Cortex=Cortex.name;
Cortex=load(Cortex);
Cortex=getfield(Cortex,xx{1});
Cortex=Cortex.*(0.195);

            if and(~contains(G{i},'trial5'),~contains(G{i},'Trial5')) %Whenever it is not PostTrial 5 

                % Sleep Scoring data
if length(states)<45*60
    states=[states nan(1,45*60-length(states))]; %Fill with NaNs.
else
    states=states(1:45*60); %Take only 45 min.
end

                %Ephys data
if length(HPC)<45*60*1000
    HPC=[HPC.' (nan(45*60*1000-length(HPC),1).')]; %Fill with NaNs.
else
    HPC=HPC(1:45*60*1000).'; %Take only 45 min.
end

if length(Cortex)<45*60*1000
    Cortex=[Cortex.' (nan(45*60*1000-length(Cortex),1).')]; %Fill with NaNs.
else
    Cortex=Cortex(1:45*60*1000).'; %Takeripple_timestamp only 45 min.
end
%      Cortex=Cortex(1:45*60*1000*4).';

%Find SD values
[sd_swr]=find_std(HPC,Cortex,states,ss);

Sd_Swr.sd2_hpc_co(i)=sd_swr.sd2_hpc_co;
Sd_Swr.sd5_hpc_co(i)=sd_swr.sd5_hpc_co;
Sd_Swr.sd2_pfc_co(i)=sd_swr.sd2_pfc_co;
Sd_Swr.sd5_pfc_co(i)=sd_swr.sd5_pfc_co;
Sd_Swr.sd2_hpc_long(i)=sd_swr.sd2_hpc_long;
Sd_Swr.sd5_hpc_long(i)=sd_swr.sd5_hpc_long;
Sd_Swr.sd2_pfc_long(i)=sd_swr.sd2_pfc_long;
Sd_Swr.sd5_pfc_long(i)=sd_swr.sd5_pfc_long;

            elseif contains(G{i}, 'rial5') % PostTrial 5 case 
                        
%Sleep scoring data
if length(states)<45*60*4
    states=[states nan(1,45*60*4-length(states))]; %Fill with NaNs.
else
    states=states(1:45*60*4); %Take only 45 min.
end


%Ephys
if length(HPC)<45*60*1000*4
    HPC=[HPC.' (nan(45*60*1000*4-length(HPC),1).')]; %Fill with NaNs.
else
    HPC=HPC(1:45*60*1000*4).'; %Take only 45 min.
end

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

HPC_1=HPC(1:2700*1000);
HPC_2=HPC(2700*1000+1:2700*2*1000);
HPC_3=HPC(1+2700*2*1000:2700*3*1000);
HPC_4=HPC(1+2700*3*1000:2700*4*1000);

Cortex_1=Cortex(1:2700*1000);
Cortex_2=Cortex(2700*1000+1:2700*2*1000);
Cortex_3=Cortex(1+2700*2*1000:2700*3*1000);
Cortex_4=Cortex(1+2700*3*1000:2700*4*1000);
%Find SD values
[sd_swr]=find_std(HPC_1,Cortex_1,states1,ss);
Sd_Swr.sd2_hpc_co(6)=sd_swr.sd2_hpc_co;
Sd_Swr.sd5_hpc_co(6)=sd_swr.sd5_hpc_co;
Sd_Swr.sd2_pfc_co(6)=sd_swr.sd2_pfc_co;
Sd_Swr.sd5_pfc_co(6)=sd_swr.sd5_pfc_co;
Sd_Swr.sd2_hpc_long(6)=sd_swr.sd2_hpc_long;
Sd_Swr.sd5_hpc_long(6)=sd_swr.sd5_hpc_long;
Sd_Swr.sd2_pfc_long(6)=sd_swr.sd2_pfc_long;
Sd_Swr.sd5_pfc_long(6)=sd_swr.sd5_pfc_long;


[sd_swr]=find_std(HPC_2,Cortex_2,states2,ss);
Sd_Swr.sd2_hpc_co(7)=sd_swr.sd2_hpc_co;
Sd_Swr.sd5_hpc_co(7)=sd_swr.sd5_hpc_co;
Sd_Swr.sd2_pfc_co(7)=sd_swr.sd2_pfc_co;
Sd_Swr.sd5_pfc_co(7)=sd_swr.sd5_pfc_co;
Sd_Swr.sd2_hpc_long(7)=sd_swr.sd2_hpc_long;
Sd_Swr.sd5_hpc_long(7)=sd_swr.sd5_hpc_long;
Sd_Swr.sd2_pfc_long(7)=sd_swr.sd2_pfc_long;
Sd_Swr.sd5_pfc_long(7)=sd_swr.sd5_pfc_long;

[sd_swr]=find_std(HPC_3,Cortex_3,states3,ss);
Sd_Swr.sd2_hpc_co(8)=sd_swr.sd2_hpc_co;
Sd_Swr.sd5_hpc_co(8)=sd_swr.sd5_hpc_co;
Sd_Swr.sd2_pfc_co(8)=sd_swr.sd2_pfc_co;
Sd_Swr.sd5_pfc_co(8)=sd_swr.sd5_pfc_co;
Sd_Swr.sd2_hpc_long(8)=sd_swr.sd2_hpc_long;
Sd_Swr.sd5_hpc_long(8)=sd_swr.sd5_hpc_long;
Sd_Swr.sd2_pfc_long(8)=sd_swr.sd2_pfc_long;
Sd_Swr.sd5_pfc_long(8)=sd_swr.sd5_pfc_long;       

[sd_swr]=find_std(HPC_4,Cortex_4,states4,ss);
Sd_Swr.sd2_hpc_co(9)=sd_swr.sd2_hpc_co;
Sd_Swr.sd5_hpc_co(9)=sd_swr.sd5_hpc_co;
Sd_Swr.sd2_pfc_co(9)=sd_swr.sd2_pfc_co;
Sd_Swr.sd5_pfc_co(9)=sd_swr.sd5_pfc_co;
Sd_Swr.sd2_hpc_long(9)=sd_swr.sd2_hpc_long;
Sd_Swr.sd5_hpc_long(9)=sd_swr.sd5_hpc_long;
Sd_Swr.sd2_pfc_long(9)=sd_swr.sd2_pfc_long;
Sd_Swr.sd5_pfc_long(9)=sd_swr.sd5_pfc_long;                                      

            end

            cd ..
        else
            cd .. %Means there is no sleep scoring file.
            end
        else
            cd ..
        end

    end
        cd ..

end


%% Threshold Table
    TT=table;
    TT.Variables=...
        [
                [{g{j}};{'x'};{'x'};{'x'};{'x'};{'x'};{'x'};{'x'}] ...
        [{'HPC_2SD_Concatenated'};{'HPC_2SD_Longest'};{'HPC_5SD_Concatenated'};{'HPC_5SD_Longest'};{'PFC_2SD_Concatenated'};{'PFC_2SD_Longest'};{'PFC_5SD_Concatenated'};{'PFC_5SD_Longest'}] ...
    num2cell([ Sd_Swr.sd2_hpc_co;Sd_Swr.sd2_hpc_long;Sd_Swr.sd5_hpc_co;Sd_Swr.sd5_hpc_long ...
        ;Sd_Swr.sd2_pfc_co;Sd_Swr.sd2_pfc_long;Sd_Swr.sd5_pfc_co;Sd_Swr.sd5_pfc_long ...
      ]) ...
    num2cell([ mean(Sd_Swr.sd2_hpc_co); mean(Sd_Swr.sd2_hpc_long); mean(Sd_Swr.sd5_hpc_co);mean(Sd_Swr.sd5_hpc_long) ...
        ;mean(Sd_Swr.sd2_pfc_co);mean(Sd_Swr.sd2_pfc_long);mean(Sd_Swr.sd5_pfc_co);mean(Sd_Swr.sd5_pfc_long) ...
      ])...
      ];
corrected_means = []; 
 TT2=cell2mat(TT{:,3:end-1});
 for m=1:size(TT2,1)
    corrected_means = [corrected_means; mean(rmoutliers(nonzeros(TT2(m,3:end))))];
 end
 TT.corrected_means = num2cell(corrected_means);
%% Iterating through trials 
       cd(g{j})
   
ripple_phase_total_c1 = {};
ripple_phase_total_c2 = {}; 
ripple_phase_total_c3 = {};
ripple_phase_total_c4 = {};

trial = [{'ps'},{'pt1'},{'pt2'},{'pt3'},{'pt4'},{'pt5.1'},{'pt5.2'},{'pt5.3'},{'pt5.4'}];

       for i = 1:length(G) % Trial Index

        cd(G{i})
        clear states
        clear HPC Cortex
        A = dir('*states*.mat');
        A = {A.name};

        if sum(contains(A, 'states')) > 0 %More than 2 sleep scoring files
            A=A(cellfun(@(x) ~isempty(strfind(x,'states')),A));
            A=A(~(cellfun(@(x) ~isempty(strfind(x,'eeg')),A)));
            if sum(contains(A, 'states')) > 0

%                     st2=st(cellfun(@(x) ~isempty(strfind(x,barea)),st)); %Brain area.
            cellfun(@load,A);



            HPC=dir('*HPC_*.mat');
            HPC=HPC.name;
            HPC=load(HPC);
            HPC=HPC.HPC;
            HPC=HPC.*(0.195);

            Cortex=dir(strcat('*',xx{1},'*.mat'));
            Cortex=Cortex.name;
            Cortex=load(Cortex);
            % Cortex=Cortex.Cortex;
            Cortex=getfield(Cortex,xx{1});
            Cortex=Cortex.*(0.195);





            if and(~contains(G{i},'trial5'),~contains(G{i},'Trial5')) %Whenever it is not PostTrial 5 

%% Phase Detection 

[swr_hpc,swr_pfc,s_hpc,s_pfc,V_hpc,V_pfc,signal2_hpc,signal2_pfc,sd_swr, M_multiplets, Mx_multiplets, multiplets, ripples, Mono_hpc, Mono_pfc,total_swrs,total_NREM_min,v_index,v_values,v_timestamps] = ripple_detection(HPC,Cortex,states,ss,offset1,offset2,TT,j,i);
% Mono_pfc is filtered

idx3 = find(contains(clus1(:,8),trial{i}));
clus1_t = clus1(idx3,:);

idx4 = find(contains(clus2(:,8),trial{i}));
clus2_t = clus2(idx4,:);

idx5 = find(contains(clus3(:,8),trial{i}));
clus3_t = clus3(idx5,:);

idx6 = find(contains(clus4(:,8),trial{i}));
clus4_t = clus4(idx6,:);

Mx_hpc_c1 = [];
Mx_hpc_c2 = [];
Mx_hpc_c3 = [];
Mx_hpc_c4 = [];
mc1 = [];
mc2 = [];
mc3 = [];
mc4 = [];

if iscell(Mono_hpc)

Mx_hpc_c1 = clus1_t(:,2);
Mx_hpc_c1 = vertcat( Mx_hpc_c1{:});
Mx_hpc_c2 = clus2_t(:,2);
Mx_hpc_c2 = vertcat(Mx_hpc_c2{:});
Mx_hpc_c3 = clus3_t(:,2);
Mx_hpc_c3 = vertcat( Mx_hpc_c3{:});
Mx_hpc_c4 = clus4_t(:,2);
Mx_hpc_c4 = vertcat( Mx_hpc_c4{:});
signal2_pfc = cellfun(@(equis) times((1/0.195), equis),Mono_pfc,'UniformOutput',false); %Remove convertion factor for ripple detection
signal2_pfc = cellfun(@(equis) mod(rad2deg(angle(hilbert(equis))),360) ,signal2_pfc,'UniformOutput',false);
ti = v_timestamps;
for ep = 1:length(v_index)
mc1{ep} = Mx_hpc_c1(Mx_hpc_c1 >= v_index(ep) & Mx_hpc_c1 <= v_index(ep)+v_values(ep)-1);
mc2{ep} = Mx_hpc_c2(Mx_hpc_c2 >= v_index(ep) & Mx_hpc_c2 <= v_index(ep)+v_values(ep)-1);
mc3{ep} = Mx_hpc_c3(Mx_hpc_c3 >= v_index(ep) & Mx_hpc_c3 <= v_index(ep)+v_values(ep)-1);
mc4{ep} = Mx_hpc_c4(Mx_hpc_c4 >= v_index(ep) & Mx_hpc_c4 <= v_index(ep)+v_values(ep)-1);
end
Mx_hpc_c1 = mc1';
Mx_hpc_c2 = mc2'; 
Mx_hpc_c3 = mc3'; 
Mx_hpc_c4 = mc4'; 

ripple_sample_c1 = cellfun(@(equis,equis2) find(ismember(equis,equis2)),ti, Mx_hpc_c1,'UniformOutput',false);
ripple_phase_c1 = cellfun(@(equis,equis2) equis(equis2),signal2_pfc,ripple_sample_c1,'UniformOutput',false);

ripple_sample_c2 = cellfun(@(equis,equis2) find(ismember(equis,equis2)),ti, Mx_hpc_c2,'UniformOutput',false);
ripple_phase_c2 = cellfun(@(equis,equis2) equis(equis2),signal2_pfc,ripple_sample_c2,'UniformOutput',false);

ripple_sample_c3 = cellfun(@(equis,equis2) find(ismember(equis,equis2)),ti, Mx_hpc_c3,'UniformOutput',false);
ripple_phase_c3 = cellfun(@(equis,equis2) equis(equis2),signal2_pfc,ripple_sample_c3,'UniformOutput',false);

ripple_sample_c4 = cellfun(@(equis,equis2) find(ismember(equis,equis2)),ti, Mx_hpc_c4,'UniformOutput',false);
ripple_phase_c4 = cellfun(@(equis,equis2) equis(equis2),signal2_pfc,ripple_sample_c4,'UniformOutput',false);
% ripple_phase_allswr = vertcat(swr_pfc{:});
% 
% ripple_phase_allclusters = [ripple_phase_c1;ripple_phase_c2];
% ripple_phase_allclusters = vertcat(ripple_phase_allclusters{:});
% 
% ripple_phase_outliers_idx = find(~ismember(ripple_phase_allswr,ripple_phase_allclusters));
% ripple_phase_outliers = ripple_phase_allswr(ripple_phase_outliers_idx);

ripple_phase_total_c1{i} = ripple_phase_c1;
ripple_phase_total_c2{i} = ripple_phase_c2;
ripple_phase_total_c3{i} = ripple_phase_c3;
ripple_phase_total_c4{i} = ripple_phase_c4;

else 
ripple_phase_total_c1{i} = NaN;
ripple_phase_total_c2{i} = NaN;
ripple_phase_total_c3{i} = NaN;
ripple_phase_total_c4{i} = NaN;

end 
            
            end

if contains(G{i}, 'rial5') %Pt5 case 
   clear HPC_1 HPC_2 HPC_3 HPC_4 Cortex_1 Cortex_2 Cortex_3 Cortex_4

                                                      
    if length(states)<45*60*4
        states=[states nan(1,45*60*4-length(states))]; %Fill with NaNs.
    else
        states=states(1:45*60*4); %Take only 45 min.
    end


    %Ephys
    if length(HPC)<45*60*1000*4
        HPC=[HPC.' (nan(45*60*1000*4-length(HPC),1).')]; %Fill with NaNs.
    else
        HPC=HPC(1:45*60*1000*4).'; %Take only 45 min.
    end

    if length(Cortex)<45*60*1000*4
        Cortex=[Cortex.' (nan(45*60*1000*4-length(Cortex),1).')]; %Fill with NaNs.
    else
    end

       for jj = 1:4
       pfc = Cortex((2700*1000*(jj-1))+1:2700*1000*jj); 
       hpc = HPC((2700*1000*(jj-1))+1:2700*1000*jj);
       states_chunk = states(2700*(jj-1)+1:2700*jj);
       %% Phase Detection for pt5 
      [swr_hpc, swr_pfc, s_hpc, s_pfc, V_hpc, V_pfc, signal2_hpc, signal2_pfc, sd_swr, M_multiplets, Mx_multiplets, multiplets, ripples, Mono_hpc, Mono_pfc,total_swrs,total_NREM_min,v_index,v_values,v_timestamps] = ripple_detection(hpc,pfc,states_chunk,ss,offset1,offset2,TT,j,i); 
        idx3 = find(contains(clus1(:,8),trial{i+jj-1}));
        clus1_t = clus1(idx3,:);

        idx4 = find(contains(clus2(:,8),trial{i+jj-1}));
        clus2_t = clus2(idx4,:);
        
        idx5 = find(contains(clus3(:,8),trial{i+jj-1}));
        clus3_t = clus3(idx5,:);
        
        idx6 = find(contains(clus4(:,8),trial{i+jj-1}));
        clus4_t = clus4(idx6,:);
        
        Mx_hpc_c1 = [];
        Mx_hpc_c2 = [];
        Mx_hpc_c3 = [];
        Mx_hpc_c4 = [];
        mc1 = [];
        mc2 = [];
        mc3 = [];
        mc4 = [];
            if iscell(Mono_hpc)

            Mx_hpc_c1 = clus1_t(:,2);
            Mx_hpc_c1 = vertcat( Mx_hpc_c1{:});
            Mx_hpc_c2 = clus2_t(:,2);
            Mx_hpc_c2 = vertcat(Mx_hpc_c2{:});
            Mx_hpc_c3 = clus3_t(:,2);
            Mx_hpc_c3 = vertcat( Mx_hpc_c3{:});
            Mx_hpc_c4 = clus4_t(:,2);
            Mx_hpc_c4 = vertcat( Mx_hpc_c4{:});
            
            signal2_pfc = cellfun(@(equis) times((1/0.195), equis),Mono_pfc,'UniformOutput',false); %Remove convertion factor for ripple detection
            signal2_pfc = cellfun(@(equis) mod(rad2deg(angle(hilbert(equis))),360) ,signal2_pfc,'UniformOutput',false);
            ti = v_timestamps;
            
            for ep = 1:length(v_index)
                mc1{ep} = Mx_hpc_c1(Mx_hpc_c1 >= v_index(ep) & Mx_hpc_c1 <= v_index(ep)+v_values(ep)-1);
                mc2{ep} = Mx_hpc_c2(Mx_hpc_c2 >= v_index(ep) & Mx_hpc_c2 <= v_index(ep)+v_values(ep)-1);
                mc3{ep} = Mx_hpc_c3(Mx_hpc_c3 >= v_index(ep) & Mx_hpc_c3 <= v_index(ep)+v_values(ep)-1);
                mc4{ep} = Mx_hpc_c4(Mx_hpc_c4 >= v_index(ep) & Mx_hpc_c4 <= v_index(ep)+v_values(ep)-1);
             end
                Mx_hpc_c1 = mc1';
                Mx_hpc_c2 = mc2'; 
                Mx_hpc_c3 = mc3'; 
                Mx_hpc_c4 = mc4'; 

            ripple_sample_c1 = cellfun(@(equis,equis2) find(ismember(equis,equis2)),ti, Mx_hpc_c1,'UniformOutput',false);
            ripple_phase_c1 = cellfun(@(equis,equis2) equis(equis2),signal2_pfc,ripple_sample_c1,'UniformOutput',false);

            ripple_sample_c2 = cellfun(@(equis,equis2) find(ismember(equis,equis2)),ti, Mx_hpc_c2,'UniformOutput',false);
            ripple_phase_c2 = cellfun(@(equis,equis2) equis(equis2),signal2_pfc,ripple_sample_c2,'UniformOutput',false);
            
            ripple_sample_c3 = cellfun(@(equis,equis2) find(ismember(equis,equis2)),ti, Mx_hpc_c3,'UniformOutput',false);
            ripple_phase_c3 = cellfun(@(equis,equis2) equis(equis2),signal2_pfc,ripple_sample_c3,'UniformOutput',false);

            ripple_sample_c4 = cellfun(@(equis,equis2) find(ismember(equis,equis2)),ti, Mx_hpc_c4,'UniformOutput',false);
            ripple_phase_c4 = cellfun(@(equis,equis2) equis(equis2),signal2_pfc,ripple_sample_c4,'UniformOutput',false);

%             ripple_phase_total_c1{i+jj-1} = ripple_phase_c1;
%             ripple_phase_total_c2{i+jj-1} = ripple_phase_c2;
%             
%             ripple_phase_allswr = vertcat(swr_pfc{:});
% 
%             ripple_phase_allclusters = [ripple_phase_c1;ripple_phase_c2];
%             ripple_phase_allclusters = vertcat(ripple_phase_allclusters{:});
% 
%             ripple_phase_outliers_idx = find(~ismember(ripple_phase_allswr,ripple_phase_allclusters));
%             ripple_phase_outliers = ripple_phase_allswr(ripple_phase_outliers_idx);

            ripple_phase_total_c1{i+jj-1} = ripple_phase_c1;
            ripple_phase_total_c2{i+jj-1} = ripple_phase_c2;
            ripple_phase_total_c3{i+jj-1} = ripple_phase_c3;
            ripple_phase_total_c4{i+jj-1} = ripple_phase_c4;

            else 
                
            ripple_phase_total_c1{i+jj-1} = NaN;
            ripple_phase_total_c2{i+jj-1} = NaN;
            ripple_phase_total_c3{i+jj-1} = NaN;
            ripple_phase_total_c4{i+jj-1} = NaN;
            
            end               

       end                                 
end
           end

        end

        cd .. 
       end
        cd ..
%% Compiling Rat wise 
ripple_phase_comp_c1 = [ripple_phase_comp_c1;ripple_phase_total_c1];
ripple_phase_comp_c2 = [ripple_phase_comp_c2;ripple_phase_total_c2];
ripple_phase_comp_c3 = [ripple_phase_comp_c3;ripple_phase_total_c3];
ripple_phase_comp_c4 = [ripple_phase_comp_c4;ripple_phase_total_c4];

%% Saving the Variables

 save(strcat('phase_clus1_4_',g{j},'.mat'),'ripple_phase_total_c1','ripple_phase_total_c2','ripple_phase_total_c3','ripple_phase_total_c4')
 save(strcat('phase_clus1_4_compilation_Rat',rat_folder{k},'.mat'),'ripple_phase_comp_c1','ripple_phase_comp_c2','ripple_phase_comp_c3','ripple_phase_comp_c4')

    end

    cd ..
        
end

%% For making individual files out of the Variables pelin sent 
 
% var_name = whos;
% var_name = {var_name.name};
% idx = find(contains(var_name, 'rat7', 'IgnoreCase',true));
% for i = 1:length(idx)
%     save(var_name{idx(i)},var_name{idx(i)})
% end
%% Plotting
%% RGS
% polarhistogram(deg2rad(ripple_phase_comp_c1_rgs(~isnan(ripple_phase_comp_c1_rgs))))
% hold on 
% polarhistogram(deg2rad(ripple_phase_comp_c2_rgs(~isnan(ripple_phase_comp_c2_rgs))))
% hold on 
% polarhistogram(deg2rad(ripple_phase_comp_c3_rgs(~isnan(ripple_phase_comp_c3_rgs))))
% hold on 
% polarhistogram(deg2rad(ripple_phase_comp_c4_rgs(~isnan(ripple_phase_comp_c4_rgs))))
% 
% legend ('Cluster1','Cluster2','Cluster3','Cluster4')
% title('Cluster wise phase analysis SWR-SO  RGS')
% 
% figure()
% 
% histogram(ripple_phase_comp_c1_rgs(~isnan(ripple_phase_comp_c1_rgs)),'Normalization','probability')
% hold on 
% histogram(ripple_phase_comp_c2_rgs(~isnan(ripple_phase_comp_c2_rgs)),'Normalization','probability')
% hold on 
% histogram(ripple_phase_comp_c3_rgs(~isnan(ripple_phase_comp_c3_rgs)),'Normalization','probability')
% hold on 
% histogram(ripple_phase_comp_c4_rgs(~isnan(ripple_phase_comp_c4_rgs)),'Normalization','probability')
% ylabel('Probability')
% legend ('Cluster1','Cluster2','Cluster3','Cluster4')
% title('Cluster wise phase analysis SWR-SO  RGS')
 %% Veh
% figure()
% polarhistogram(deg2rad(ripple_phase_comp_c1_veh(~isnan(ripple_phase_comp_c1_veh))),'EdgeColor','b', 'LineWidth',1)
% hold on 
% polarhistogram(deg2rad(ripple_phase_comp_c2_veh(~isnan(ripple_phase_comp_c2_veh))),'EdgeColor','k', 'LineWidth',1)
% hold on 
% polarhistogram(deg2rad(ripple_phase_comp_c3_veh(~isnan(ripple_phase_comp_c3_veh))),'EdgeColor','g', 'LineWidth',1)
% hold on 
% polarhistogram(deg2rad(ripple_phase_comp_c4_veh(~isnan(ripple_phase_comp_c4_veh))),'EdgeColor','r', 'LineWidth',1)
% 
% legend ('Cluster1','Cluster2','Cluster3','Cluster4')
% title('Cluster wise phase analysis SWR-SO  Veh')
% 
% figure()
% 
% histogram(ripple_phase_comp_c1_veh(~isnan(ripple_phase_comp_c1_veh)),'Normalization','probability','BinWidth',4,'EdgeColor','b', 'LineWidth',1.5)
% hold on 
% histogram(ripple_phase_comp_c2_veh(~isnan(ripple_phase_comp_c2_veh)),'Normalization','probability','BinWidth',4,'EdgeColor','k', 'LineWidth',1.5)
% hold on 
% histogram(ripple_phase_comp_c3_veh(~isnan(ripple_phase_comp_c3_veh)),'Normalization','probability','BinWidth',4,'EdgeColor','g', 'LineWidth',1.5)
% hold on 
% histogram(ripple_phase_comp_c4_veh(~isnan(ripple_phase_comp_c4_veh)),'Normalization','probability', 'BinWidth',4,'EdgeAlpha',0.2,'EdgeColor','r', 'LineWidth',1.5)
% ylabel('Probability')
% legend ('Cluster1','Cluster2','Cluster3','Cluster4')
% title('Cluster wise phase analysis SWR-SO  Veh')