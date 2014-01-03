function [points,normals] = cloud_process(points,normals,H)
k=6;
lambda=10;
n=size(points,1);
dim=size(points,2);

tic
[i,j,v,E]=compute_AE(points,normals,k,3); % k变化了，需要重新计算
E=E+1;
A=sparse(i+1,j+1,v,n*dim*k,n*dim);
t=toc;
fprintf(1,'compute AE takes:%f\n',t);

% ==========================================================================================
% ADMM solver, no ||N_in^i-N_out^i||\leq gamma constraints

% normals = normalOpt(points, normals, lambda, k, false, A , E, H);
% normals = normalize_normals(normals);

% reweighted
for i=1:1
    normals = normalOpt(points, normals, lambda, k, true, A, E, H);
    normals = normalize_normals(normals);
end

% ===========================================================================================
write_mesh(points,normals,'out_normal.xyzn');
return;

% position optimization
lambda=10;
[points] = optimizePos(points, normals, k, lambda, A, E);

% output
write_mesh(points,normals,'out.xyzn');
