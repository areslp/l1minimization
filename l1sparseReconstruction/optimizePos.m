function [Xout] = optimizePos(Xin, Nout, k, lambda)

N = size( Xin, 1 );
dim=size(Nout,2);

%% 使用k-nearest方法求得E
E = zeros( k * N, 2 ); % E存储Nin中索引
tree = kdtree_build( Xin );
for i = 1:1:N
    idxs = kdtree_k_nearest_neighbors( tree, Xin( i, : ), k+1 );
    idxs(idxs==i)=[];
    E( k*(i-1)+1:1:k*i, 1 ) = i;
    E( k*(i-1)+1:1:k*i, 2 ) = idxs;
end
fprintf('finish calculating E %dx%d\n', size(E,1), size(E,2));

%% 计算W
% sigma = 10 / 180 * pi;
% W = zeros( k*N, 1 );
% for i = 1:1:k*N
    % p1 = Nout( E( i, 1), :);
    % p2 = Nout( E( i, 2), :);
    % theta = acos( dot(p1, p2) / ( norm(p1, 2) * norm(p2, 2) ) );
    % W( i ) = exp( - (theta / sigma)^4 );
% end
% fprintf('finish calculating W %dx%d\n', size(W,1), size(W,2));

%% 计算N_avg,E中两点向量均值
Navg = zeros(k * N, dim);
for i = 1:1:k*N
    n1 = Nout(E(i, 1), :);
    n2 = Nout(E(i, 2), :);
    Navg(i, :) = (n1 + n2) ./ 2;
end
fprintf('finish calculating Navg %dx%d\n', size(Navg,1), size(Navg,2));

%% 计算A
% A = sparse( [1:1:k*N, 1:1:k*N], [E(:,1); E(:,2)],...
% [W.*(sum(Navg.*Nout(E(:,1), :), 2)), -W.*(sum(Navg.*Nout(E(:,2), :), 2))], k*N, N);
% fprintf('finish calculating A %dx%d\n', size(A,1), size(A,2));

%% 计算f
% f = sum(Navg .* (Xin(E(:,1), :) - Xin(E(:,2), :)), 2);
% fprintf('finish calculating f %dx%d\n', size(f,1), size(f,2));

% min ||At+f||_1+lambda||t||_2
% t=l1minimization(A,f,lambda); % lambda越大，变化越小

D=kdtree_adj_vec(Xin,k,dim);
b=reshape(Xin',dim*N,1); % vectorize
x=l1minimization2(Navg,D,b,lambda); % lambda越大，变化越小

%% 得出Xout
% Xout = Xin + repmat(t, 1, 3) .* Nout;
Xout=reshape(x,dim,N)';
fprintf('finish calculating Xout %dx%d\n', size(Xout,1), size(Xout,2));
end
