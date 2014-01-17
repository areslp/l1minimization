clear all;
close all;

[points,normals]=mesh_import('paper_data/noise_0.5_fandisk2.xyzn');

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
