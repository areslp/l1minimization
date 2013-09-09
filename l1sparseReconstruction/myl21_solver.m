function [ Nout,x ] = l1SparseReconstruct( points, normals, k, lambda )
n=size( points, 1 );
N=reshape(normals',3*n,1); % vectorize
B=construct_Adj(points,normals,k,false); % B is the adjacent matrix
%% 最后的方程求解
% minimize 1/2*|| Ax - b ||_2^2 + \lambda sum(norm(x_i))

A=speye(3*k*n);
b=B*N;

% verify B and B1, correct
% tmp=reshape(b,3,k*n)';
% norm(B1*Nin-tmp)

norm(B*N,1)

[x] = group_lasso_wrapper(B,A,b,lambda);
% size(x)
norm(x,1)

[NO]=lsqr_wrapper(B,x);

norm(B*NO,1)

% norm(NO-N)

Nout=reshape(NO,3,n)';

error=zeros(n,1);
for i=1:n
    error(i)=norm(normals(i,:)-Nout(i,:));
end
max(error)
