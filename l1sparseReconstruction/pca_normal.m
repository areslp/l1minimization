function [normals] = pca_normal(points,k)
n=size(points,1);
dim=size(points,2);
normals=zeros(n,dim);
tree = kdtree_build( points);
for i=1:n
    idxs = kdtree_k_nearest_neighbors( tree, points( i, : ), k+1 );
    X=points(idxs,:);
    COEFF = pca(X);
    normals(i,:)=COEFF(:,2);
end
