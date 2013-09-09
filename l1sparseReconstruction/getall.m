clear all;
close all;

data=load('fandisk.xyzn');
% data=load('noise_0.5_fandisk.xyzn');
% data=load('noise_0.5_fandisk2.xyzn');
% data=load('noise_outliers_0.5_fandisk.xyzn');
points=data(:,1:3);
normals=data(:,4:6);
normals = normals./repmat(sqrt(sum(normals.^2,2)),1,size(normals,2));
k=6;
% lambda=0.25; % 无噪声
lambda=0.25; % 0.5噪声


% correct normal
% cdata=load('standard_normal.xyzn');
% cpoints=cdata(:,1:3);
% cnormals=cdata(:,4:6);
% cnormals = cnormals./repmat(sqrt(sum(cnormals.^2,2)),1,size(cnormals,2));
% disp('correct data');
% B=construct_Adj_cvx(cpoints,cnormals,k,false); % B is the adjacent matrix
% S=B*cnormals;
% obj=sqrt(sum(S.^2,2));
% figure;
% hist(obj,20);
% title('correct normal');
% obj=sum(obj);
% fprintf(2,'obj of correct normal is %f\n',obj);


% disp('input data');
% B=construct_Adj_cvx(points,normals,k,false); % B is the adjacent matrix
% S=B*normals;
% obj=sqrt(sum(S.^2,2));
% figure;
% hist(obj,20);
% title('input normal');
% obj=sum(obj);
% fprintf(2,'obj of pca normal is %f\n',obj);

% fprintf(1,'point size is %d\n',size(points,1));

% before optimization
% verify_result(points,normals,normals,k);

% norm(normals-cnormals)

% ==========================================================================================

[normals_new,x]=myl21_solver(points,normals,k,lambda);
% disp('norm(cx-x)');
% norm(cx-x)
% normalize
normals_new = normals_new./repmat(sqrt(sum(normals_new.^2,2)),1,size(normals_new,2));
% verify_result
% verify_result(points,normals,normals_new,k);

% second iteration, weighted
% lambda=0.07;
normals=normals_new;
normals_new=myl21_solver_weighted(points,normals,k,lambda);
% normalize
normals_new = normals_new./repmat(sqrt(sum(normals_new.^2,2)),1,size(normals_new,2));
% verify_result
% verify_result(points,normals,normals_new,k);
% TODO: 结果有可能有Inf，查查哪里的问题
inan=isnan(normals_new);
normals_new(inan(:,1),:)=[];
points(inan(:,1),:)=[];

iinf=isinf(normals_new);
normals_new(iinf(:,1),:)=[];
points(iinf(:,1),:)=[];

save NormalStep.mat normals_new points;

% ===========================================================================================
ff1=fopen('out_normal.xyzn','w');
for i=1:length(points)
    fprintf(ff1,'%f %f %f %f %f %f\n',points(i,1),points(i,2),points(i,3),normals_new(i,1),normals_new(i,2),normals_new(i,3));
end
fclose(ff1);

return;

% cvx version

% TODO: using the correct normal for testing

% points=cpoints;
% normals=cnormals;

% first iteration
% disp('first iteration');
% % lambda=0.13;
% lambda=0.8;
% normals_new=cvx_wrapper(points,normals,k,lambda,true); % cvx version
% normals_new = normals_new./repmat(sqrt(sum(normals_new.^2,2)),1,size(normals_new,2));
% normals=normals_new;

% B=construct_Adj_cvx(points,normals,k,false); % B is the adjacent matrix
% S=B*normals;
% obj=sum(S.^2,2);
% figure;
% hist(obj,20);
% title('first iteration, optimized normal');
% obj=sum(obj);
% fprintf(2,'obj of first iteration is %f\n',obj);

% second iteration
% disp('second iteration');
% lambda=0.05;
% normals_new=cvx_wrapper(points,normals,k,lambda,true); % cvx version
% normals_new = normals_new./repmat(sqrt(sum(normals_new.^2,2)),1,size(normals_new,2));
% normals=normals_new;

% B=construct_Adj_cvx(points,normals,k,true); % B is the adjacent matrix
% S=B*normals;
% obj=sum(S.^2,2);
% figure;
% hist(obj,20);
% title('second iteration, optimized normal');
% obj=sum(obj);
% fprintf(2,'obj of second iteration is %f\n',obj);

normals=normals_new;
normals=full(normals);

% position optimization
[points] = optimizePos(points, normals, k);

% output
ff=fopen('out.xyzn','w');
for i=1:length(points)
    fprintf(ff,'%f %f %f %f %f %f\n',points(i,1),points(i,2),points(i,3),normals(i,1),normals(i,2),normals(i,3));
end
fclose(ff);
