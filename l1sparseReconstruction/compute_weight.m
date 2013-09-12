function [W] = compute_weight(points,normals,E)
n=size(points,1);
dim=size(points,2);
k=size(E,1)/(n*dim);
% compute W
sigma = 10 / 180 * pi;
W = zeros(dim*k*n,1);
index=1; % 处理第几条边
for i = 1:dim:dim*k*n
    idx1=E(i,1);
    idx2=E(i,2);
    if rem(idx1,dim)==0
        idx1=floor(idx1/dim)-1;
    else
        idx1=floor(idx1/dim);
    end
    idx1=idx1+1;
    if rem(idx2,dim)==0
        idx2=floor(idx2/dim)-1;
    else
        idx2=floor(idx2/dim);
    end
    idx2=idx2+1;
    n1 = normals(idx1,:);
    n2 = normals(idx2,:);
    theta = acos( dot(n1, n2) / ( norm(n1, 2) * norm(n2, 2) ) );
    weight=exp(-(theta/sigma)^4);
    W(dim*(index-1)+1:dim*(index-1)+dim) = weight;
    index=index+1;
end
W=sparse(W);
