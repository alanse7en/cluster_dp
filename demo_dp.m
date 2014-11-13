%% demo of clustering by finding density peaks

clear all
close all
clc
% Load data and calculate distance matrix
data = load('fig2_panelB.dat');
dist = pdist2(data,data);
%  Clustering
para.method = 'gaussian';
para.percent = 2.0;
[cluster_lables, center_idxs] = cluster_dp(dist, para);
% show the result
figure(1);
cmap = colormap;
nclust = length(center_idxs);
for i = 1:nclust
    tmp_data = data(cluster_lables==i,:);
    ic = int8((i*64.)/(nclust*1.));
    col = cmap(ic,:);
    plot(tmp_data(:,1),tmp_data(:,2),...
        'o','MarkerSize',2,'MarkerFaceColor',col,'MarkerEdgeColor',col);
    hold on;
end
tmp_data = data(cluster_lables==0,:);
plot(tmp_data(:,1),tmp_data(:,2),...
        'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');