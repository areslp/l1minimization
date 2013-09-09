function [ Nout ] = l1SparseReconstruct( Xin, Nin, k )

n = size( Xin, 1 );

%% 使用k-nearest方法求得E
tic;
E = zeros(3*k*n,2); % E存储Nin中索引
tree = kdtree_build( Xin );
for i = 1:n
    idxs = kdtree_k_nearest_neighbors( tree, Xin( i, : ), k );
    for j=1:k
        idx=idxs(j); % the idx point
        % ith point start from 3*k*(i-1)

        % the g_index for i
        g_index_i=3*(i-1)+1;
        g_index_idx=3*(idx-1)+1;
        % the g_index for idx
        E(3*k*(i-1)+3*(j-1)+1,1)=g_index_i;
        E(3*k*(i-1)+3*(j-1)+2,1)=g_index_i+1;
        E(3*k*(i-1)+3*(j-1)+3,1)=g_index_i+2;
        E(3*k*(i-1)+3*(j-1)+1,2)=g_index_idx;
        E(3*k*(i-1)+3*(j-1)+2,2)=g_index_idx+1;
        E(3*k*(i-1)+3*(j-1)+3,2)=g_index_idx+2;
    end
end
kdtree_delete(tree);
toc;
N=reshape(Nin',3*n,1);

%% 计算W
tic;
sigma = 10 / 180 * pi;
W = zeros(3*k*n,1);
for i = 1:k*n
    idx1=E(i,1);
    idx2=E(i,2);

    if rem(idx1,3)==0
        idx1=floor(idx1/3)-1;
    else
        idx1=floor(idx1/3);
    end
    idx1=idx1+1;

    if rem(idx2,3)==0
        idx2=floor(idx2/3)-1;
    else
        idx2=floor(idx2/3);
    end
    idx2=idx2+1;

    p1 = Nin(idx1,:);
    p2 = Nin(idx2,:);
    theta = acos( dot(p1, p2) / ( norm(p1, 2) * norm(p2, 2) ) );
    weight=exp(-(theta/sigma)^4);
    W(3*(i-1)+1) = weight;
    W(3*(i-1)+2) = weight;
    W(3*(i-1)+3) = weight;
end
W = sparse(W);
toc;


tic;
size(E)
3*k*n
A = sparse( [1:1:3*k*n, 1:1:3*k*n], [E(:,1); E(:,2)],...
    [ones(3*k*n, 1), -ones(3*k*n, 1)], 3*k*n, n);
toc;

%% 最后的方程求解

% 1/2||ZM-N'||_F^2+lambda||X||_2,1+mu/2||X-Z+Y/mu||_F^2
lambda=0.1;
b=W*A*N;
rho=1.1;
alpha=1.0;
p=ones(n,1);
p=p*3;
[z, history] = group_lasso(speye(3*n), b, lambda, p, rho, alpha);
Nout=reshape(z,3,n)';
