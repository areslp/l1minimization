function [Nout] = normalOpt(points, normals, lambda, k, weighted, A, E, H);
n = size( points, 1 );
dim = size(normals,2);
N=reshape(normals',dim*n,1); % vectorize
D=A;
if weighted
    W=compute_weight(points,normals,E);
    D=diag(W)*A;
end
%% 最后的方程求解
% minimize 1/2*|| x - b ||_2^2 + \lambda ||Dx||_1

% v0=compute_Dx(D,N,3);
% [x]=total_variation_vec(N,lambda,3,D);
% [x]=total_variation_vec_H(N,lambda,3,D,H);
% [x]=tvl1_total_variation_vec_H(N,lambda,D,H);
[x]=tvl1_total_variation_vec(N,lambda,D);
% [x]=tvl2_total_variation_vec(N,lambda,D);
% [x]=tvl2_total_variation_vec_H(N,lambda,D,H);
% v=compute_Dx(D,x,3);
% fprintf(2,'init Dx is %f\n',v0);
% fprintf(2,'L1 Dx is %f\n',v);
Nout=reshape(x,dim,n)';
