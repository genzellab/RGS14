%% histograms - ripple type - raw/bp - HC/OS
clc
clear
% load data of condition (MANUALLY)
addpath(genpath('/Users/pelinozsezer/Documents/MATLAB/Toolboxes')); % needs padcat


%% extract features
data=waveforms_cluster1_raw_veh; % MANUALLY CHANGE!

% HC/OS split
HC_i = contains(data(:,8),'HC');
OS_i = ~contains(data(:,8),'HC');
HC = data(HC_i,:);
OS = data(OS_i,:);

% features
si=[];
for i=1:length(HC(:,1))
    duration_start = 3001-(HC{i,3}-HC{i,2})*1000;
    duration_end   = 3001+(HC{i,4}-HC{i,3})*1000;
    si{i}=(HC{i}(2,floor(duration_start):round(duration_end)));
end
si=si';

si2=[];
for i=1:length(OS(:,1))
    duration_start = 3001-(OS{i,3}-OS{i,2})*1000;
    duration_end   = 3001+(OS{i,4}-OS{i,3})*1000;
    si2{i}=(OS{i}(2,floor(duration_start):round(duration_end)));
end
si2=si2';

timeasleep=0;
print_hist=1;

[x,y,z,w,h,q,l,p,si_mixed,th,PCA_features]=delta_specs(si,timeasleep,print_hist); % input=waveform
[x2,y2,z2,w2,h2,q2,l2,p2,si_mixed2,th2,PCA_features2]=delta_specs(si2,timeasleep,print_hist);

HC=PCA_features(:,2:end);  % check PCA_features!!
OS=PCA_features2(:,2:end); 


%% 1. mean freq (Hz)
% enter the data
data1=HC(:,1);
data2=OS(:,1);

% PLOTTING 
data= padcat(data1,data2);
% p-value
[h,p] = kstest2(data1,data2);
close all

h2= histogram(data1,'Normalization','probability');
h2.NumBins=125;
h2.FaceAlpha=1;
h2.FaceColor='r';
hold on

h1=histogram(data2,'Normalization','probability');
h1.NumBins=125;
h1.FaceAlpha=0.5;
h1.FaceColor='b';

formatSpec = '%.4f';
str_handle=["Veh: mean freq - c1 - HC vs OS (p = " + num2str(p,formatSpec) + ")"];
title(str_handle);
ylabel('Probability');
legend('HC','OS');
xlabel('mean freq');

saveas(gcf,'meanfreq_c1_raw_veh.jpg')
saveas(gcf,'meanfreq_c1_raw_veh.pdf')
close all


%% 2. amplitude (mV)
% enter the data
data1=HC(:,2);
data2=OS(:,2);

% PLOTTING 
data= padcat(data1,data2);
% p-value
[h,p] = kstest2(data1,data2);
close all

h2= histogram(data1,'Normalization','probability');
h2.NumBins=125;
h2.FaceAlpha=1;
h2.FaceColor='r';
hold on

h1=histogram(data2,'Normalization','probability');
h1.NumBins=125;
h1.FaceAlpha=0.5;
h1.FaceColor='b';
hold on

formatSpec = '%.4f';
str_handle=["Veh: amplitude - c1 - HC vs OS (p = " + num2str(p,formatSpec) + ")"];
title(str_handle);
ylabel('Probability');
legend('HC','OS');
xlabel('Amplitude');

saveas(gcf,'amplitude_c1_raw_veh.jpg')
saveas(gcf,'amplitude_c1_raw_veh.pdf')
close all

%% 3. area under the curve
% enter the data
data1=HC(:,3);
data2=OS(:,3);

% PLOTTING 
data= padcat(data1,data2);
% p-value
[h,p] = kstest2(data1,data2);
close all

h2= histogram(data1,'Normalization','probability');
h2.NumBins=125;
h2.FaceAlpha=1;
h2.FaceColor='r';
hold on

h1=histogram(data2,'Normalization','probability');
h1.NumBins=125;
h1.FaceAlpha=0.5;
h1.FaceColor='b';
hold on


formatSpec = '%.4f';
str_handle=["Veh: area und. the curve - c1 - HC vs OS (p = " + num2str(p,formatSpec) + ")"];
title(str_handle);
ylabel('Probability');
legend('HC','OS');
xlabel('Area under the curve');

saveas(gcf,'areaunder_c1_raw_veh.jpg')
saveas(gcf,'areaunder_c1_raw_veh.pdf')
close all


%% 4. duration (s/ms)
% enter the data
data1=HC(:,4);
data2=OS(:,4);

% PLOTTING 
data= padcat(data1,data2);
% p-value
[h,p] = kstest2(data1,data2);
close all

h2= histogram(data1,'Normalization','probability');
h2.NumBins=20;
h2.FaceAlpha=1;
h2.FaceColor='r';
hold on

h1=histogram(data2,'Normalization','probability');
h1.NumBins=20;
h1.FaceAlpha=0.5;
h1.FaceColor='b';
hold on

formatSpec = '%.4f';
str_handle=["Veh: duration - c1 - HC vs OS (p = " + num2str(p,formatSpec) + ")"];
title(str_handle);
ylabel('Probability');
legend('HC','OS');
xlabel('Duration');

saveas(gcf,'duration_c1_raw_veh.jpg')
saveas(gcf,'duration_c1_raw_veh.pdf')
close all


%% 5. peak 2 peak distance (mV)
% enter the data
data1=HC(:,5);
data2=OS(:,5);

% PLOTTING 
data= padcat(data1,data2);
% p-value
[h,p] = kstest2(data1,data2);
close all

h2= histogram(data1,'Normalization','probability');
h2.NumBins=125;
h2.FaceAlpha=1;
h2.FaceColor='r';
hold on

h1=histogram(data2,'Normalization','probability');
h1.NumBins=125;
h1.FaceAlpha=0.5;
h1.FaceColor='b';
hold on

formatSpec = '%.4f';
str_handle=["Veh: peak2peak dist. - c1 - HC vs OS (p = " + num2str(p,formatSpec) + ")"];
title(str_handle);
ylabel('Probability');
legend('HC','OS');
xlabel('peak2peak distance');

saveas(gcf,'peak2peak_c1_raw_veh.jpg')
saveas(gcf,'peak2peak_c1_raw_veh.pdf')
close all


%% 6. power
% enter the data
data1=HC(:,6);
data2=OS(:,6);

% PLOTTING 
data= padcat(data1,data2);
% p-value
[h,p] = kstest2(data1,data2);
close all

h2= histogram(data1,'Normalization','probability');
h2.NumBins=125;
h2.FaceAlpha=1;
h2.FaceColor='r';
hold on

h1=histogram(data2,'Normalization','probability');
h1.NumBins=125;
h1.FaceAlpha=0.5;
h1.FaceColor='b';
hold on

formatSpec = '%.4f';
str_handle=["Veh: power - c1 - HC vs OS (p = " + num2str(p,formatSpec) + ")"];
title(str_handle);
ylabel('Probability');
legend('HC','OS');
xlabel('Power');

saveas(gcf,'power_c1_raw_veh.jpg')
saveas(gcf,'power_c1_raw_veh.pdf')
close all


%% 7. entropy
% enter the data
data1=HC(:,7);
data2=OS(:,7);

% PLOTTING 
data= padcat(data1,data2);
% p-value
[h,p] = kstest2(data1,data2);
close all

h2= histogram(data1,'Normalization','probability');
h2.NumBins=10;
h2.FaceAlpha=1;
h2.FaceColor='r';
hold on


h1=histogram(data2,'Normalization','probability');
h1.NumBins=10;
h1.FaceAlpha=0.5;
h1.FaceColor='b';
hold on

formatSpec = '%.4f';
str_handle=["Veh: entropy - c1 - HC vs OS (p = " + num2str(p,formatSpec) + ")"];
title(str_handle);
ylabel('Probability');
legend('HC','OS');
xlabel('Entropy');

saveas(gcf,'entropy_c1_raw_veh.jpg')
saveas(gcf,'entropy_c1_raw_veh.pdf')
close all


%% 8. number of peaks
% enter the data
data1=HC(:,8);
data2=OS(:,8);

% PLOTTING 
data= padcat(data1,data2);
% p-value
[h,p] = kstest2(data1,data2);
close all

h2= histogram(data1,'Normalization','probability');
h2.NumBins=20;
h2.FaceAlpha=1;
h2.FaceColor='r';
hold on

h1=histogram(data2,'Normalization','probability');
h1.NumBins=20;
h1.FaceAlpha=0.5;
h1.FaceColor='b';
hold on

formatSpec = '%.4f';
str_handle=["Veh: # of peaks - c1 - HC vs OS (p = " + num2str(p,formatSpec) + ")"];
title(str_handle);
ylabel('Probability');
legend('HC','OS');
xlabel('# of peaks');

saveas(gcf,'noofpeaks_c1_raw_veh.jpg')
saveas(gcf,'noofpeaks_c1_raw_veh.pdf')
close all




