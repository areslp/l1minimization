clear all;
close all;

[points,normals]=mesh_import('noise_0.5_fandisk2.xyzn');
% [points,normals]=mesh_import('noise_0.5_mechpart.xyzn');
% [points,normals]=mesh_import('noise_0.5_oilpmp.xyzn');
% [points,normals]=mesh_import('noise_0.5_screwdriver.xyzn'); % 螺丝刀没法处理，太尖
% [points,normals]=mesh_import('noise_0.1_maxplanck.xyzn');
% [points,normals]=mesh_import('noise_0.2_bunny.xyzn');
% [points,normals]=mesh_import('noise_0.5_double-torus2.xyzn');
% [points,normals]=mesh_import('ToyTurtle.xyzn');
% [points,normals]=mesh_import('noise_0.1_hand.xyzn');

%% smooth normal

% n=size(points,1);
% dim=size(points,2);
% Asize=n*dim;
% smooth_kernel=30;
% tic
% [i,j,v]=compute_WL(points,smooth_kernel);
% H=sparse(i+1,j+1,v,Asize,Asize);
% t=toc;
% fprintf(1,'compute blur kernel takes:%f\n',t);

% vectorize normals
% normal_vector=reshape(normals',Asize,1); % vectorize
% normal_vector=H*normal_vector;
% normals=reshape(normal_vector,dim,n)';
% normals = normalize_normals(normals);
% write_mesh(points,normals,'blurred_mesh.xyzn');
% [points,normals]=mesh_import('blurred_mesh.xyzn');

% cloud_process(points,normals,H);
cloud_process(points,normals,zeros(1));
