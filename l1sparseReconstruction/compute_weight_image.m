function [D] = compute_weight(I,E)
rowe=size(E,1);
[m n k]=size(I);
dim=k;
W = zeros(rowe,1);
V=reshape(I,m*n,k);
sigma=0.5;
index=1;
for i = 1:dim:rowe
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
    v1 = V(idx1,:);
    v2 = V(idx2,:);
    weight=exp(-(norm(v1-v2)/sigma)^4);
    W(dim*(index-1)+1:dim*(index-1)+dim) = weight;
    index=index+1;
end
W=sparse(W);
B = sparse( [1:rowe, 1:rowe], [E(:,1); E(:,2)],...
    [ones(rowe, 1), -ones(rowe, 1)], rowe, m*n*k);
D=diag(W)*B;
