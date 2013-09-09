function [ Nout ] = l1SparseReconstruct( Xin, Nin, k )

n = size( Xin, 1 );
Nout = zeros(n, 3);

%% 使用k-nearest方法求得E
tic;
E = zeros( k * n, 2 ); % E存储Nin中索引
tree = kdtree_build( Xin );
for i = 1:1:n
    idxs = kdtree_k_nearest_neighbors( tree, Xin( i, : ), k );
    E( k*(i-1)+1:1:k*i, 1 ) = i;
    E( k*(i-1)+1:1:k*i, 2 ) = idxs;
end
kdtree_delete(tree);
toc;
% test
% figure;
% hold on;
% t = 10;
% T = randperm( N, t );
% for i = 1:1:t
%     p = (T(i)-1)*k+1;
%     a = E( p, 1 );
%     scatter3( Xin(a, 1), Xin(a, 2), Xin(a, 3), 'r.' );
%     b = E(p:1:p+k-1, 2);
%     scatter3( Xin(b, 1), Xin(b, 2), Xin(b, 3), 'b.' );
% end

%% 计算W
tic;
sigma = 10 / 180 * pi;
W = zeros( k*n, 1 );
for i = 1:1:k*n
    p1 = Nin( E( i, 1), :);
    p2 = Nin( E( i, 2), :);
    theta = acos( dot(p1, p2) / ( norm(p1, 2) * norm(p2, 2) ) );
    W( i ) = exp( - (theta / sigma)^4 );
end
W = sparse(W);
toc;
%Nout = W;

%% 求A，(ni - nj) = A * Nin
%{
A = zeros( k * N, N );
for i = 1:1:k*N
    A( i, E(i, 1) ) = 1;
    A( i, E(i, 2) ) = -1;
end
%}
tic;
A = sparse( [1:1:k*n, 1:1:k*n], [E(:,1); E(:,2)],...
    [ones(k*n, 1), -ones(k*n, 1)], k*n, n);
toc;

% test
% T1 = A * Nin;
% T2 = Nin( E(:,1), : ) - Nin(E(:,2), :);
% res = T1 - T2
%% 最后的方程求解

% 1/2||ZM-N'||_F^2+lambda||X||_2,1+mu/2||X-Z+Y/mu||_F^2
tic
B = diag(W) * A; % kn*n
[u s v] = svds(B);
Binv = v / s * u'; % n*kn
t=toc;
fprintf(1,'the inv of B takes %f\n',t);

% precompute vars
tic;
M=Binv'; % kn*n
mt=M'; % n*kn
nt=Nin'; % 3*n
mmt=M*M'; % kn*kn
ntmt=nt*mt; % 3*kn
kn=k*n;
t=toc;
fprintf(1,'precomputed vars takes %f\n',t);

% optimization vars
X=zeros(3,k*n);
Z=X;

% parameters
Y=X;
lambda = 0.1;
mu = 1e-6;
iter=1;
max_iter=1000;
while true 
    tic
    cfQ=Z-Y/mu;
    cflambda=lambda/mu;
    X = l21(cfQ, cflambda);
    t=toc;
    fprintf(1,'L21 takes %f\n',t);
    tic
    Z = (ntmt/mu+X+Y/mu)/(mmt/mu+speye(kn));
    t=toc;
    fprintf(1,'update Z(inv) takes %f\n',t);
    tic;
    Y = Y + mu * (X - Z);
    mu = min(1.1 * mu, 1e10);
    iter=iter+1;
    t1=norm(X-Z,1);
    t2=norm(X*M-nt,1);
    if iter==1 || mod(iter,50)==0 
        fprintf(1,'iter %d,reconstruction error is %f,change error is %f\n',iter,t2,t1);
    end
    if t1 <= 10^-6 && t2 <= 10^-6
        break;
        fprintf(1,'converged, iter %d,reconstruction error is %f,change error is %f\n',iter,t2,t1);
    end
    if iter>max_iter
        break;
        fprintf(1,'max iteration reached! t1 is %f,t2 is %f\n',t1,t2);
    end
    t=toc;
    fprintf(1,'others takes %f\n',t);
end;

Nout = M * X';

