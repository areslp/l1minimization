function [B] = construct_A(Xin,Nin,k,weighted) 
n=size(Xin,1);
dim=size(Nin,2);
%% 使用k-nearest方法求得E
% tic;
E = zeros(dim*k*n,2); % E存储Nin中索引
tree = kdtree_build( Xin );
for i = 1:n
    idxs = kdtree_k_nearest_neighbors( tree, Xin( i, : ), k );
    for j=1:k
        idx=idxs(j); % the idx point
        % ith point start from 3*k*(i-1)
        % the g_index for i
        g_index_i=dim*(i-1)+1;
        g_index_idx=dim*(idx-1)+1;
        % the g_index for idx
        for tt=1:dim
            E(dim*k*(i-1)+dim*(j-1)+tt,1)=g_index_i+tt-1;
            E(dim*k*(i-1)+dim*(j-1)+tt,2)=g_index_idx+tt-1;
        end
    end
end
kdtree_delete(tree);
% toc;

if weighted
    %% 计算W
    % tic;
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
        n1 = Nin(idx1,:);
        n2 = Nin(idx2,:);
        theta = acos( dot(n1, n2) / ( norm(n1, 2) * norm(n2, 2) ) );
        weight=exp(-(theta/sigma)^4);
        for tt=1:dim
            W(dim*(index-1)+tt) = weight;
        end
        index=index+1;
    end
    W=sparse(W);
end

A = sparse( [1:1:dim*k*n, 1:1:dim*k*n], [E(:,1); E(:,2)],...
    [ones(dim*k*n, 1), -ones(dim*k*n, 1)], dim*k*n, dim*n);

if weighted
    B=diag(W)*A; 
else
    B=A;
end
