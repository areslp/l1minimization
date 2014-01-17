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

% [x]=total_variation_vec(N,lambda,3,D);
[x]=tvl1_total_variation_vec(N,lambda,D); %% 大论文使用
% [x]=tvl2_total_variation_vec(N,lambda,D);
Nout=reshape(x,dim,n)';
