function [x] = Landweber_iteration(A,b,AtA,iter)
% min 1/2||Ax-b||_2^2
Atb=A'*b;
[m,n]=size(A);
x=zeros(n,1);
for i=1:iter
    x=x-AtA*x+Atb;
end
