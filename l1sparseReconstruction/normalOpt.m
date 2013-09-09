function [Nout] = normalOpt(points, normals, lambda, k, weighted, varargin);
n = size( points, 1 );
dim = size(normals,2);
N=reshape(normals',dim*n,1); % vectorize
if weighted
    D=construct_Adj(points,normals,k,true);
else
    D=construct_Adj(points,normals,k,false); % B is the adjacent matrix
end
%% 最后的方程求解
% minimize 1/2*|| x - b ||_2^2 + \lambda ||Dx||_1

norm(D*N,1)

% [x]=total_variation(N,lambda,D);
[x]=total_variation_vec(N,lambda,dim,D,k);

norm(D*x,1)

Nout=reshape(x,dim,n)';
