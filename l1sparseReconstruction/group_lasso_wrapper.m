function [x] = group_lasso_wrapper(B,A,b,lambda)
rho=1.1;
alpha=1.0;
[m,n]=size(B);
kn=m/3;
p=ones(kn,1);
p=p*3;
[x, history] = group_lasso(A, b, lambda, p, rho, alpha);

% obj1=1/2*norm(A*x-b)*norm(A*x-b);
% q=reshape(x,3,kn)';
% sparsev=sqrt(sum(q.^2,2));
% obj2=norm(sparsev,1);
% fprintf(1,'obj1 is %f,obj2 is %f,total obj is %f\n',obj1,obj2,obj1+lambda*obj2);
% fprintf(1,'sparsity is %f\n',nnz(sparsev)/length(sparsev));
