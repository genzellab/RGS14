
%% Pelin Özsezer
% Extracts ripples from waveforms;
% Generates the data of features;
% Normalizes and computes PCA;
% Generates density maps; Thresholds;
% Fits GMM model; Finds clusters per treatment;
% Creates new files per cluster and per treatment.

clc
clear
format compact
format longG
load('GC_window_treatmentwise.mat'); % load data

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EXTRACT WAVEFORMS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%% Vehicle
GC_bandpassed_veh = GC_window_ripples_comp_veh(:,:);
GC_raw_veh = GC_window_ripples_broadband_comp_veh(:,:);
I=cellfun(@isempty,GC_raw_veh(:,1));
I=not(I);
GC_raw_veh=GC_window_ripples_broadband_comp_veh(I,:);
GC_bandpassed_veh = GC_window_ripples_comp_veh(I,:);

L=cellfun(@length, GC_raw_veh(:,1),'UniformOutput', false);
L=cell2mat(L);
GC_raw_veh=GC_raw_veh(L~=1,:);
GC_bandpassed_veh=GC_bandpassed_veh(L~=1,:);

data=[];
for i=1:length(GC_raw_veh(:,1))
    duration_start = 3001-(GC_raw_veh{i,3}-GC_raw_veh{i,2})*1000;
    duration_end   = 3001+(GC_raw_veh{i,4}-GC_raw_veh{i,3})*1000;
    data{i}=(GC_raw_veh{i}(2,floor(duration_start):round(duration_end)));
end
waveforms_veh=data';

%% RGS

GC_bandpassed_rgs = GC_window_ripples_comp_rgs(:,:);
GC_raw_rgs = GC_window_ripples_broadband_comp_rgs(:,:);
I=cellfun(@isempty,GC_raw_rgs(:,1));
I=not(I);
GC_raw_rgs=GC_window_ripples_broadband_comp_rgs(I,:);
GC_bandpassed_rgs = GC_window_ripples_comp_rgs(I,:);

L=cellfun(@length, GC_raw_rgs(:,1),'UniformOutput', false);
L=cell2mat(L);
GC_raw_rgs=GC_raw_rgs(L~=1,:);
GC_bandpassed_rgs=GC_bandpassed_rgs(L~=1,:);

data=[];
for i=1:length(GC_raw_rgs(:,1))
    duration_start = 3001-(GC_raw_rgs{i,3}-GC_raw_rgs{i,2})*1000;
    duration_end   = 3001+(GC_raw_rgs{i,4}-GC_raw_rgs{i,3})*1000;
    data{i}=(GC_raw_rgs{i}(2,floor(duration_start):round(duration_end)));
end
waveforms_rgs=data';

%%%%%%%%%%%%%%%%
%%% FEATURES %%%
%%%%%%%%%%%%%%%%
data=waveforms_veh;
si=data;
data2=waveforms_rgs;
si2=data2;

timeasleep=0;
print_hist=1;

[x,y,z,w,h,q,l,p,si_mixed,th,PCA_features]=delta_specs(si,timeasleep,print_hist); % Vehicle
[x2,y2,z2,w2,h2,q2,l2,p2,si_mixed2,th2,PCA_features2]=delta_specs(si2,timeasleep,print_hist); % RGS

PCA_features_veh_m22=PCA_features(:,2:end);
PCA_features_rgs_m22=PCA_features2(:,2:end);

%%%%%%%%%%%%%%%%%%%%%
%%% COMPUTING PCA %%%
%%%%%%%%%%%%%%%%%%%%%
X=zscore(PCA_features_veh_m22);
Y=zscore(PCA_features_rgs_m22);

[coeff,score,latent,~,explained] = pca(X);
[coeff2,score2,latent2,~,explained2] = pca(Y);

Z_veh=X*coeff;
Z_rgs=Y*coeff;
%Z_rgs=Y*coeff2;

%%%%%%%%%%%%%%%%%%%%%%%
%%% 3D DENSITY MAPS %%%
%%%%%%%%%%%%%%%%%%%%%%%
% histcn NEEDED!
% https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/23897/versions/4/previews/histcn.m/index.html
addpath(genpath('/Users/pelinozsezer/Documents/MATLAB/Toolboxes'));


% Select the data you want to analyze! (Veh or RGS)
% Veh
X=Z_veh(:,1); % Change the value according to which PCA/s you want to analyze!
Y=Z_veh(:,2); % e.g.; 1 corresponds to PCA1, etc.
Z=Z_veh(:,3);
% RGS
X=Z_rgs(:,1);
Y=Z_rgs(:,2);
Z=Z_rgs(:,3);

data = [X,Y,Z];

%% Needed part of the script of heatscatter()
numbins = 125;  % change numbin HERE
markersize = 10;
marker = 'o';
plot_colorbar = 1;
plot_lsf = 1;

%% values
%[values, centers] = hist3([X Y], [numbins numbins]);
[values centers mid loc] = histcn(data);
values = values/size(X,1); % to produce density

centers_X = centers{1,1};
centers_Y = centers{1,2};
centers_Z = centers{1,3};

binsize_X = abs(centers_X(2) - centers_X(1)) / 2;
binsize_Y = abs(centers_Y(2) - centers_Y(1)) / 2;
binsize_Z = abs(centers_Z(2) - centers_Z(1)) / 2;
bins_X = zeros(numbins, 2);
bins_Y = zeros(numbins, 2);
bins_Z = zeros(numbins, 2);

for i = 1:numbins
    bins_X(i, 1) = centers_X(i) - binsize_X;
    bins_X(i, 2) = centers_X(i) + binsize_X;
    bins_Y(i, 1) = centers_Y(i) - binsize_Y;
    bins_Y(i, 2) = centers_Y(i) + binsize_Y;
    bins_Z(i, 1) = centers_Z(i) - binsize_Z;
    bins_Z(i, 2) = centers_Z(i) + binsize_Z;
end

scatter_COL = zeros(length(X), 1);

onepercent = round(length(X) / 100);

fprintf('Generating colormap...\n');

for i = 1:length(X)

    if (mod(i,onepercent) == 0)
        fprintf('.');
    end

    last_lower_X = NaN;
    last_higher_X = NaN;
    id_X = NaN;

    c_X = X(i);
    last_lower_X = find(c_X >= bins_X(:,1));
    if (~isempty(last_lower_X))
        last_lower_X = last_lower_X(end);
    else
        last_higher_X = find(c_X <= bins_X(:,2));
        if (~isempty(last_higher_X))
            last_higher_X = last_higher_X(1);
        end
    end
    if (~isnan(last_lower_X))
        id_X = last_lower_X;
    else
        if (~isnan(last_higher_X))
            id_X = last_higher_X;
        end
    end

    last_lower_Y = NaN;
    last_higher_Y = NaN;
    id_Y = NaN;

    c_Y = Y(i);
    last_lower_Y = find(c_Y >= bins_Y(:,1));
    if (~isempty(last_lower_Y))
        last_lower_Y = last_lower_Y(end);
    else
        last_higher_Y = find(c_Y <= bins_Y(:,2));
        if (~isempty(last_higher_Y))
            last_higher_Y = last_higher_Y(1);
        end
    end
    if (~isnan(last_lower_Y))
        id_Y = last_lower_Y;
    else
        if (~isnan(last_higher_Y))
            id_Y = last_higher_Y;
        end
    end

    last_lower_Z = NaN;
    last_higher_Z = NaN;
    id_Z = NaN;

    c_Z = Z(i);
    last_lower_Z = find(c_Z >= bins_Z(:,1));
    if (~isempty(last_lower_Z))
        last_lower_Z = last_lower_Z(end);
    else
        last_higher_Z = find(c_Z <= bins_Z(:,2));
        if (~isempty(last_higher_Z))
            last_higher_Z = last_higher_Z(1);
        end
    end
    if (~isnan(last_lower_Z))
        id_Z = last_lower_Z;
    else
        if (~isnan(last_higher_Z))
            id_Z = last_higher_Z;
        end
    end

    scatter_COL(i) = values(id_X, id_Y, id_Z);

end

%%%%%%%%%%%%%%%%%%%%
%%% THRESHOLDING %%%
%%%%%%%%%%%%%%%%%%%%

%% THRESHOLD - change values accordingly
thresholded_idx =find(scatter_COL>0.0005);  % change threshold value here!
%
scatter_COL=scatter_COL(thresholded_idx,:);
X=X(thresholded_idx,:);
Y=Y(thresholded_idx,:);
Z=Z(thresholded_idx,:);

%     f = figure();
%     scatter3(X, Y,Z, markersize, scatter_COL, marker);
%     xlabel('PCA1');
%     ylabel('PCA2');
%     zlabel('PCA3');
%     title('')
%     hold on
%
%     if (plot_colorbar)
%         cbh = colorbar;
%     end
%     colormap("hot");
%
%     saveas(gcf,'thresholded_pca.jpg');
%     saveas(gcf,'thresholded_pca.pdf');

%% Find the best model!
AIC = zeros(1,10);
GMModels = cell(1,10);
options = statset('MaxIter',500);
for k = 1:10
    GMModels{k} = fitgmdist(data,k,'Options',options,'CovarianceType','diagonal');
    AIC(k)= GMModels{k}.AIC;
end
AIC=AIC';

BIC = zeros(1,10);
GMModels = cell(1,10);
options = statset('MaxIter',500);
for k = 1:10
    GMModels{k} = fitgmdist(data,k,'Options',options,'CovarianceType','diagonal');
    BIC(k)= GMModels{k}.BIC;
end
BIC=BIC';

plot(AIC,'-r')
hold on
plot(BIC,'-b')
legend('AIC','BIC')
xlabel('number of clusters')

%% Find the best number of component with kmeans!
rng('default') % For reproducibility
eva = evalclusters(data,'kmeans','CalinskiHarabasz','KList',1:6);
eva.OptimalK

%%%%%%%%%%%%%%%
%%% FIT GMM %%%
%%%%%%%%%%%%%%%

data = [X,Y,Z];

%% fit the model
rng(1);
k =4; % change the number of components/clusters here!
options = statset('MaxIter',1000);
gmm_data = fitgmdist(data,k,'CovarianceType','diagonal','SharedCovariance',false,'Options',options);
cluster_data=cluster(gmm_data,data);

%% cluster1
idx = cluster_data==1; % CHANGE NUMBER ACCORDINGLY!
actual = data(idx,:);

x_axis1 =actual(:,1);
y_axis1 =actual(:,2);
z_axis1 =actual(:,3);

rect_x_points1=[mean(x_axis1)-2*std(x_axis1) mean(x_axis1)+2*std(x_axis1)];
rect_y_points1=[mean(y_axis1)-2*std(y_axis1) mean(y_axis1)+2*std(y_axis1)];
rect_z_points1=[mean(z_axis1)-2*std(z_axis1) mean(z_axis1)+2*std(z_axis1)];

width1     = rect_x_points1(2)-rect_x_points1(1);
height1    = rect_y_points1(2)-rect_y_points1(1);
depth1     = rect_z_points1(2)-rect_z_points1(1);

%% cluster 2
idx = cluster_data==2; % CHANGE NUMBER ACCORDINGLY!
actual = data(idx,:);

x_axis2 =actual(:,1);
y_axis2 =actual(:,2);
z_axis2 =actual(:,3);

rect_x_points2=[mean(x_axis2)-2*std(x_axis2) mean(x_axis2)+2*std(x_axis2)];
rect_y_points2=[mean(y_axis2)-2*std(y_axis2) mean(y_axis2)+2*std(y_axis2)];
rect_z_points2=[mean(z_axis2)-2*std(z_axis2) mean(z_axis2)+2*std(z_axis2)];

width2     = rect_x_points2(2)-rect_x_points2(1);
height2    = rect_y_points2(2)-rect_y_points2(1);
depth2     = rect_z_points2(2)-rect_z_points2(1);

%% cluster 3
idx = cluster_data==3; % CHANGE NUMBER ACCORDINGLY!
actual = data(idx,:);

x_axis3 =actual(:,1);
y_axis3 =actual(:,2);
z_axis3 =actual(:,3);

rect_x_points3=[mean(x_axis3)-2*std(x_axis3) mean(x_axis3)+2*std(x_axis3)];
rect_y_points3=[mean(y_axis3)-2*std(y_axis3) mean(y_axis3)+2*std(y_axis3)];
rect_z_points3=[mean(z_axis3)-2*std(z_axis3) mean(z_axis3)+2*std(z_axis3)];

width3     = rect_x_points3(2)-rect_x_points3(1);
height3    = rect_y_points3(2)-rect_y_points3(1);
depth3     = rect_z_points3(2)-rect_z_points3(1);

%% cluster 4
idx = cluster_data==4; % CHANGE NUMBER ACCORDINGLY!
actual = data(idx,:);

x_axis4 =actual(:,1);
y_axis4 =actual(:,2);
z_axis4 =actual(:,3);

rect_x_points4=[mean(x_axis4)-2*std(x_axis4) mean(x_axis4)+2*std(x_axis4)];
rect_y_points4=[mean(y_axis4)-2*std(y_axis4) mean(y_axis4)+2*std(y_axis4)];
rect_z_points4=[mean(z_axis4)-2*std(z_axis4) mean(z_axis4)+2*std(z_axis4)];

width4     = rect_x_points4(2)-rect_x_points4(1);
height4    = rect_y_points4(2)-rect_y_points4(1);
depth4     = rect_z_points4(2)-rect_z_points4(1);


%% Plot
scatter3(x_axis1,y_axis1,z_axis1,'filled','MarkerFaceColor',[209/255 233/255 196/255],'MarkerEdgeColor',[209/255 233/255 196/255])
hold on
scatter3(x_axis2,y_axis2,z_axis2,'filled','MarkerFaceColor',[44/255 210/255 245/255],'MarkerEdgeColor',[44/255 210/255 245/255])
hold on
scatter3(x_axis3,y_axis3,z_axis3,'filled','MarkerFaceColor',[10/255 75/255 141/255],'MarkerEdgeColor',[10/255 75/255 141/255])
hold on
scatter3(x_axis4,y_axis4,z_axis4,'filled','MarkerFaceColor',[0 120/255 0],'MarkerEdgeColor',[0 120/255 0])
hold on
title('4 clusters');
xlabel('PCA1')
ylabel('PCA2')
zlabel('PCA3')
legend('cluster 1','cluster 2','cluster 3','cluster 4');
xlim([-10 60])
ylim([-30 20])
zlim([-10 40])

saveas(gcf,'4clusters.jpg')
saveas(gcf,'4clusters.pdf')
saveas(gcf,'4clusters.fig')


%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CROP CLUSTER DATA %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%% cluster 1
cluster_1_idx=[];
idx_no=[];
for i=1:length(x_axis1)
    sprintf('x_axis   %.0f', i)
    idx_no=find(X==x_axis1(i));
    cluster_1_idx=[cluster_1_idx; idx_no] ;
end

%% cluster 2
cluster_2_idx=[];
idx_no=[];
for i=1:length(x_axis2)
    sprintf('x_axis   %.0f', i)
    idx_no=find(X==x_axis2(i));
    cluster_2_idx=[cluster_2_idx; idx_no] ;
end

%% cluster 3
cluster_3_idx=[];
idx_no=[];
for i=1:length(x_axis3)
    sprintf('x_axis   %.0f', i)
    idx_no=find(X==x_axis3(i));
    cluster_3_idx=[cluster_3_idx; idx_no] ;
end

%% cluster 4
cluster_4_idx=[];
idx_no=[];
for i=1:length(x_axis4)
    sprintf('x_axis   %.0f', i)
    idx_no=find(X==x_axis4(i));
    cluster_4_idx=[cluster_4_idx; idx_no] ;
end


%% Veh
cluster1_idx_veh=cluster_1_idx;
cluster2_idx_veh=cluster_2_idx;
cluster3_idx_veh=cluster_3_idx;
cluster4_idx_veh=cluster_4_idx;
%thresholded_idx_veh=thresholded_idx;

%% RGS
cluster1_idx_rgs=cluster_1_idx;
cluster2_idx_rgs=cluster_2_idx;
cluster3_idx_rgs=cluster_3_idx;
cluster4_idx_rgs=cluster_4_idx;
%thresholded_idx_rgs=thresholded_idx;


%% Veh
waveforms_cluster1_raw_veh=GC_raw_veh;
waveforms_cluster1_raw_veh= waveforms_cluster1_raw_veh(cluster1_idx_veh,:);

waveforms_cluster2_raw_veh=GC_raw_veh;
waveforms_cluster2_raw_veh= waveforms_cluster2_raw_veh(cluster2_idx_veh,:);

waveforms_cluster3_raw_veh=GC_raw_veh;
waveforms_cluster3_raw_veh= waveforms_cluster3_raw_veh(cluster3_idx_veh,:);

waveforms_cluster4_raw_veh=GC_raw_veh;
waveforms_cluster4_raw_veh= waveforms_cluster4_raw_veh(cluster4_idx_veh,:);

%

waveforms_cluster1_bp_veh=GC_bandpassed_veh;
waveforms_cluster1_bp_veh= waveforms_cluster1_bp_veh(cluster1_idx_veh,:);

waveforms_cluster2_bp_veh=GC_bandpassed_veh;
waveforms_cluster2_bp_veh= waveforms_cluster2_bp_veh(cluster2_idx_veh,:);

waveforms_cluster3_bp_veh=GC_bandpassed_veh;
waveforms_cluster3_bp_veh= waveforms_cluster3_bp_veh(cluster3_idx_veh,:);

waveforms_cluster4_bp_veh=GC_bandpassed_veh;
waveforms_cluster4_bp_veh= waveforms_cluster4_bp_veh(cluster4_idx_veh,:);


%% RGS
waveforms_cluster1_raw_rgs=GC_raw_rgs;
waveforms_cluster1_raw_rgs= waveforms_cluster1_raw_rgs(cluster1_idx_rgs,:);

waveforms_cluster2_raw_rgs=GC_raw_rgs;
waveforms_cluster2_raw_rgs= waveforms_cluster2_raw_rgs(cluster2_idx_rgs,:);

waveforms_cluster3_raw_rgs=GC_raw_rgs;
waveforms_cluster3_raw_rgs= waveforms_cluster3_raw_rgs(cluster3_idx_rgs,:);

waveforms_cluster4_raw_rgs=GC_raw_rgs;
waveforms_cluster4_raw_rgs= waveforms_cluster4_raw_rgs(cluster4_idx_rgs,:);

%

waveforms_cluster1_bp_rgs=GC_bandpassed_rgs;
waveforms_cluster1_bp_rgs= waveforms_cluster1_bp_rgs(cluster1_idx_rgs,:);

waveforms_cluster2_bp_rgs=GC_bandpassed_rgs;
waveforms_cluster2_bp_rgs= waveforms_cluster2_bp_rgs(cluster2_idx_rgs,:);

waveforms_cluster3_bp_rgs=GC_bandpassed_rgs;
waveforms_cluster3_bp_rgs= waveforms_cluster3_bp_rgs(cluster3_idx_rgs,:);

waveforms_cluster4_bp_rgs=GC_bandpassed_rgs;
waveforms_cluster4_bp_rgs= waveforms_cluster4_bp_rgs(cluster4_idx_rgs,:);


%%
clearvars -except waveforms_cluster1_raw_veh waveforms_cluster2_raw_veh waveforms_cluster3_raw_veh waveforms_cluster4_raw_veh...
    waveforms_cluster1_bp_veh waveforms_cluster2_bp_veh waveforms_cluster3_bp_veh waveforms_cluster4_bp_veh...
    waveforms_cluster1_raw_rgs waveforms_cluster2_raw_rgs waveforms_cluster3_raw_rgs waveforms_cluster4_raw_rgs...
    waveforms_cluster1_bp_rgs waveforms_cluster2_bp_rgs waveforms_cluster3_bp_rgs waveforms_cluster4_bp_rgs

save('waveforms_ripple_clusters_all.mat','-v7.3')


