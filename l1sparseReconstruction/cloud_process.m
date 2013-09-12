function [points,normals] = cloud_process(points,normals)
k=3;
lambda=0.2;

[A,E]=kdtree_adj_vec(points,k,3); % k变化了，需要重新计算
% ==========================================================================================
% ADMM solver, no ||N_in^i-N_out^i||\leq gamma constraints
% normals = normalOpt(points, normals, lambda, k, false, A , E);
% normals = normalize_normals(normals);
lambda=10;
for i=1:1
    normals = normalOpt(points, normals, lambda, k, true, A, E);
    normals = normalize_normals(normals);
end

% ===========================================================================================
ff1=fopen('out_normal.xyzn','w');
for i=1:length(points)
    fprintf(ff1,'%f %f %f %f %f %f\n',points(i,1),points(i,2),points(i,3),normals(i,1),normals(i,2),normals(i,3));
end
fclose(ff1);

return;

% position optimization
lambda=0.2;
[points] = optimizePos(points, normals, k, lambda, A, E);

% output
ff=fopen('out.xyzn','w');
for i=1:length(points)
    fprintf(ff,'%f %f %f %f %f %f\n',points(i,1),points(i,2),points(i,3),normals(i,1),normals(i,2),normals(i,3));
end
fclose(ff);
