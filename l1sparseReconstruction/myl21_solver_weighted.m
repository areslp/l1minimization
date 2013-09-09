function [ Nout ] = l1SparseReconstruct( Xin, Nin, k, lambda )

n = size( Xin, 1 );

N=reshape(Nin',3*n,1);

B=construct_Adj(Xin,Nin,k,true);

%% 最后的方程求解
% minimize 1/2*|| Ax - b ||_2^2 + \lambda sum(norm(x_i))
A=speye(3*k*n);
b=B*N;

[x] = group_lasso_wrapper(B,A,b,lambda);

[NO]=lsqr_wrapper(B,x);

norm(NO-N)

Nout=reshape(NO,3,n)';
