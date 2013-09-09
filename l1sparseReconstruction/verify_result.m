function [x] = verify_result(points,normals,normals_new,k)

n = size(points , 1);

B=construct_Adj_cvx(points,normals_new,k,false); % B is the adjacent matrix

V=B*normals_new;
sparsev=sqrt(sum(V.^2,2));
% sparsev=sum(V,2);

% max(sparsev)
% nnz(sparsev)
% length(sparsev)

fprintf(2,'the object value is %f\n',norm(sparsev,1));
fprintf(2,'sparsity is %f\n',nnz(sparsev)/length(sparsev));
fprintf(2,'norm(n-n_new) is %f\n',norm(normals_new-normals,'fro'));

% 3*k*n
x=reshape(V',3*k*n,1); % vectorize
