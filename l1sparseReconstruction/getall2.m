clear all;
close all;

% [points,normals]=mesh_import('ori_data/double-torus1.xyzn');
[points,normals]=mesh_import('standard_normal.xyzn');
% [points,normals]=mesh_import('fandisk.xyzn');

n=size(points,1);
dim=size(points,2);
Asize=n*dim;
smooth_kernel=9;

tic
[i,j,v]=compute_WL(points,smooth_kernel);
H=sparse(i+1,j+1,v,Asize,Asize);
t=toc;
fprintf(1,'compute blur kernel takes:%f\n',t);

% vectorize normals
normal_vector=reshape(normals',Asize,1); % vectorize
normal_vector=H*normal_vector;
normals=reshape(normal_vector,dim,n)';
normals = normalize_normals(normals);
write_mesh(points,normals,'blurred_mesh.xyzn');

[points,normals]=mesh_import('blurred_mesh.xyzn');

[null_normal_idx,col]=find(isnan(normals)==1);
points(null_normal_idx,:)=[];
normals(null_normal_idx,:)=[];
write_mesh(points,normals,'not_null_mesh.xyzn');
cloud_process(points,normals,H);
% cloud_process2(points,normals,H);
