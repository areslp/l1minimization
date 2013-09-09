function [normals] = normalize_normals(normals)
normals = normals./repmat(sqrt(sum(normals.^2,2)),1,size(normals,2));
