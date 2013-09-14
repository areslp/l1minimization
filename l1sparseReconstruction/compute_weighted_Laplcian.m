function [A] = compute_weighted_Laplcian(points,k)
% 构造mesh的带权拉普拉斯矩阵
n=size(points,1);
dim=size(points,2);
Asize=n*dim;
tree = kdtree_build( points);
A=sparse([],[],[],Asize,Asize,Asize*(k+1));
index=1; % 表示当前处理到A的index行
sigma=0.1;
for i=1:n
    idxs=kdtree_k_nearest_neighbors( tree, points( i, : ), k+1 );
    p=points(i,:);
    pmat=repmat(p,k+1,1);
    dmat=points(idxs,:);
    % 计算i点与idxs点之间的距离
    dis=sparse([],[],[],n,1,k+1);
    dis(idxs)=sqrt(sum((pmat-dmat).^2,2));
    % 生成权重向量
    dis(idxs)=exp(-dis(idxs).^2/sigma^2);
    dis(idxs)=dis(idxs)/norm(dis(idxs));
    pidx=sparse([],[],[],n,1,k+1);
    pidx(idxs)=1;
    psf=pidx.*dis;
    % 生成A矩阵的三行
    psf=repmat(psf',3,1);
    psf=reshape(psf,1,3*n);
    A(i,:)=psf;
end
kdtree_delete(tree);
