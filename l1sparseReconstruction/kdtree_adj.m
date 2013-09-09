function [A] = kdtree_adj(points,k)

n = size( points, 1 );

tree = kdtree_build( points );
%% 使用k-nearest方法求得E
% tic;
E = zeros( n*k, 2 ); % E存储Nin中索引
index=1;
for i = 1:1:n
    idxs = kdtree_k_nearest_neighbors( tree, points( i, : ), k+1 );
    idxs(find(idxs==i))=[];
    nn=length(idxs);
    E( index:index+nn-1, 1 ) = i;
    E( index:index+nn-1, 2 ) = idxs;
    index=index+nn;
end
kdtree_delete(tree);

enum=length(E); % edge number

A = sparse( [1:1:enum, 1:1:enum], [E(:,1); E(:,2)],...
    [ones(enum, 1), -ones(enum, 1)], enum, n);
