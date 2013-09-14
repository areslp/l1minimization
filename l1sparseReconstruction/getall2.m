clear all;
close all;

% data=load('double-torus2.xyzn');
% data=load('consist_noise_0.5_double-torus1.xyzn');
% data=load('noise_0.5_fandisk.xyzn');
data=load('ori_data/double-torus1.xyzn');


points=data(:,1:3);
normals=data(:,4:6);
normals = normals./repmat(sqrt(sum(normals.^2,2)),1,size(normals,2));

tic
H=compute_weighted_Laplcian(points,normals,6);
toc
pause;

cloud_process(points,normals);
