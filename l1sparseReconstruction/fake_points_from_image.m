function [points] = fake_points_from_image(I)
[m,n,k]=size(I);
% 第一列
a=[1:m];
a=repmat(a,n,1);
a=reshape(a,m*n,1);
b=[1:n]';
b=repmat(b,m,1);
c=zeros(m*n,1);
points=[a b c];
