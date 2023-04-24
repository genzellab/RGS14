cd('/home/adrian/Documents/GitHub/RGS14/spiking_activity')

% RGS14
%load('Name_Confirming_Data.mat');
load('New_Wake_Data_Extra_Units_clean.mat');

x=New_Wake_Data_Extra_Units;
%% 
v=[];
v_name={};
for i=1:length(x) % 107 pyr units
    
    v(i,:)=[x(i).Presleep x(i).Trial_1 x(i).Post_Trial_1 x(i).Trial_2 x(i).Post_Trial_2 x(i).Trial_3 x(i).Post_Trial_3 x(i).Trial_4 x(i).Post_Trial_4 x(i).Trial_5 x(i).Post_Trial_5_All ];
    v_name{i}=x(i).WFM_Title;
end
% %% Getting correct dates. 
% y=Name_Confirming_Data;
% v_confirm={};
% for i=1:length(y)
%         v_confirm{i}=y(i).WFM_Titles;
% end

%% Veh units
load('Veh_Final_Units_clean.mat')
x=Veh_Final_Units;
%%
v=[];
v_name={};
for i=1:length(x) % 57 pyr units
    
    v(i,:)=[x(i).Presleep x(i).Trial_1 x(i).Post_Trial_1 x(i).Trial_2 x(i).Post_Trial_2 x(i).Trial_3 x(i).Post_Trial_3 x(i).Trial_4 x(i).Post_Trial_4 x(i).Trial_5 x(i).Post_Trial_5_All ];
    v_name{i}=x(i).WFM_Title;
end

