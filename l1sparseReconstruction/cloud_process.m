function [points,normals] = cloud_process(points,normals)
k=6;
% lambda=0; % 无噪声
lambda=0.2; % 0.5噪声

% ==========================================================================================
% ADMM solver, no ||N_in^i-N_out^i||\leq gamma constraints
normals = normalOpt(points, normals, lambda, k, false);
normals = normalize_normals(normals);
lambda=0.1;
for i=1:2
    normals = normalOpt(points, normals, lambda, k, true);
    normals = normalize_normals(normals);
end

% CVX solver, the same as paper
% [ normals_new] = cvx_wrapper(points, normals, k, lambda, false);
% normals_new = normalize_normals(normals_new);
% normals=normals_new;
% [normals_new] = cvx_wrapper(points, normals, k, lambda, true);
% normals_new = normalize_normals(normals_new);

% ===========================================================================================
ff1=fopen('out_normal.xyzn','w');
for i=1:length(points)
    fprintf(ff1,'%f %f %f %f %f %f\n',points(i,1),points(i,2),points(i,3),normals(i,1),normals(i,2),normals(i,3));
end
fclose(ff1);

% save NormalStep.mat normals points;
% load NormalStep.mat normals points;
return;

% position optimization
k=6;
lambda=0.1;
[points] = optimizePos(points, normals, k, lambda);

% output
ff=fopen('out.xyzn','w');
for i=1:length(points)
    fprintf(ff,'%f %f %f %f %f %f\n',points(i,1),points(i,2),points(i,3),normals(i,1),normals(i,2),normals(i,3));
end
fclose(ff);
