function [B] = construct_A(Xin,Nin,k,weighted) 
n = size( Xin, 1 );
Nout = zeros(n, 3);
tree = kdtree_build( Xin );
%% 使用k-nearest方法求得E
% tic;
E=zeros(k*n,2); % E存储Nin中索引
index=1;
for i=1:n
    idxs = kdtree_k_nearest_neighbors( tree, Xin( i, : ), k+1 );
    idxs(find(idxs==i))=[];
    nn=length(idxs);
    E(index:index+nn-1, 1) = i;
    E(index:index+nn-1, 2) = idxs;
    index=index+nn;
end
kdtree_delete(tree);
enum=length(E); % edge number
%% 计算W
% tic;
sigma = 10 / 180 * pi;
W = zeros( enum, 1 );
for i = 1:1:enum
    p1 = Nin( E( i, 1), :);
    p2 = Nin( E( i, 2), :);
    theta = acos( dot(p1, p2) / ( norm(p1, 2) * norm(p2, 2) ) );
    W( i ) = exp( - (theta / sigma)^4 );
end
W = sparse(W);
% toc;
%Nout = W;

%% 求A，(ni - nj) = A * Nin
%{
A = zeros( k * N, N );
for i = 1:1:k*N
    A( i, E(i, 1) ) = 1;
    A( i, E(i, 2) ) = -1;
end
%}
% tic;
A = sparse( [1:1:enum, 1:1:enum], [E(:,1); E(:,2)],...
    [ones(enum, 1), -ones(enum, 1)], enum, n);
% toc;

if weighted
    B=diag(W)*A;
else
    B=A;
end
