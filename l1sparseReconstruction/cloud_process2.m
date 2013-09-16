function [] = cloud_process2(points,normals,k)
% 先解z=Dx, min 1/2||z-DN||_2^2+lambda||z||_1
% 然后再解min 1/2||z-Dx||_2^2+lambda||x-y||_2^2

k=3;
lambda=0.5;
n=size(points,1);
dim=size(points,2);

tic
[i,j,v,E]=compute_AE(points,normals,k,3); % k变化了，需要重新计算
E=E+1;
A=sparse(i+1,j+1,v,n*dim*k,n*dim);
t=toc;
fprintf(1,'compute AE takes:%f\n',t);

% ==========================================================================================
% ADMM solver
y=reshape(normals',dim*n,1); % vectorize
D=A;

% 解z
Dy=D*y;
A=speye(dim*k*n);
[z] = group_lasso_wrapper(D,A,Dy,lambda);
% 解x，z=Dx min 1/2||x-y||_2^2 s.t. z=Dx
% lambda=1e-5;
% x=constraint_least_square(D,y,z,lambda);
[x]=lsqr_wrapper(D,z);

Nout=reshape(x,dim,n)';

% ===========================================================================================
write_mesh(points,normals,'out_normal.xyzn');
return;

