function [H] = compute_WL_1D(points,k)
% 类似compute_WL，但是计算的是对点云里一维使用的模糊核
n=size(points,1);
%% 使用k-nearest方法求得E
tree = kdtree_build( points);
H=sparse([],[],[],n,n,n*(k+1));
sigma=1;
for i = 1:n
    idxs = kdtree_k_nearest_neighbors( tree, points( i, : ), k+1 );
    sum_v=0;
    for j=1:(k+1)
        jidx=idxs(j);
        dis=norm(points(i,:)-points(jidx,:));
        v=exp(-dis*dis/(sigma*sigma));
        sum_v=sum_v+v;
    end
    for j=1:(k+1)
        jidx=idxs(j);
        dis=norm(points(i,:)-points(jidx,:));
        v=exp(-dis*dis/(sigma*sigma));
        H(i,jidx)=v/sum_v;
    end
end
kdtree_delete(tree);
