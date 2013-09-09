function [A] = kdtree_adj(points,k,vec)
% 对三维数据用kdtree构造邻接矩阵
% assert(vec==3);
dim=vec;
n=size(points,1);
%% 使用k-nearest方法求得E
% tic;
E = zeros(dim*k*n,2); % E存储Nin中索引
tree = kdtree_build( points);
for i = 1:n
    idxs = kdtree_k_nearest_neighbors( tree, points( i, : ), k+1 );
    idxs(idxs==i)=[];
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

A = sparse( [1:1:dim*k*n, 1:1:dim*k*n], [E(:,1); E(:,2)],...
    [ones(dim*k*n, 1), -ones(dim*k*n, 1)], dim*k*n, dim*n);
