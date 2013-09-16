function [points,normals] = mesh_import(file)
data=load(file);
points=data(:,1:3);
normals=data(:,4:6);
normals = normals./repmat(sqrt(sum(normals.^2,2)),1,size(normals,2));
