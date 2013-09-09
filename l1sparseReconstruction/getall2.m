clear all;
close all;

% data=load('double-torus1.xyzn');
data=load('double-torus2.xyzn');
% data=load('consist_out.xyzn');
% data=load('noise_0.2_double-torus2.xyzn');
% data=load('fandisk.xyzn');
% data=load('noise_0.5_fandisk.xyzn');
% data=load('consist_noise_0.5_fandisk.xyzn'); % 当PCA邻域取很小时，有可能法矢方向不一致
% data=load('noise_0.5_fandisk2.xyzn');
% data=load('noise_outliers_0.5_fandisk.xyzn');
% data=load('standard_normal.xyzn');


points=data(:,1:3);
normals=data(:,4:6);
normals = normals./repmat(sqrt(sum(normals.^2,2)),1,size(normals,2));

cloud_process(points,normals);
