clear all;
close all;

k=4;
lambda=0.2;

[points,normals]=gen_test();
normals = normalize_normals(normals);
n=size(points,1);
dim=size(points,2);
smooth_kernel=9;
ori_normals=normals;
Asize=n*dim;
[A,E]=kdtree_adj_vec(points,k,3); % k变化了，需要重新计算
% normals=pca_normal(points,k); % TODO: 不是consist的法矢
% draw_points_and_normals(points,normals,'ori');
N=reshape(normals',dim*n,1); % vectorize
vc=compute_Dx(A,N,3);
fprintf(2,'correct Dx: %f\n',vc);

tic
[i,j,v]=compute_WL(points,smooth_kernel);
H=sparse(i+1,j+1,v,Asize,Asize);
t=toc;
fprintf(1,'compute blur kernel takes:%f\n',t);

% vectorize normals
normal_vector=reshape(normals',Asize,1); % vectorize
normal_vector=H*normal_vector;
normals=reshape(normal_vector,dim,n)';
normals = normalize_normals(normals);
draw_points_and_normals(points,normals,'blurred');

% 分维度plot原始信号和blur后的信号
% for i=1:dim
    % b0=ori_normals(:,i);
    % b=normals(:,i);
    % figure;
    % hold on;
    % plot( b, 'r-o' );
    % plot( b0, 'b-x' );
    % hold off;
% end
% pause;

% 将normal分解为3个1维的子信号，分别在上面求解tv
% A = sparse( -diag( ones(n,1), 0 ) + diag( ones(n-1,1), 1 ) );
% A(n,n-1)=1; % use a backward difference for the final point
% H1D=compute_WL_1D(points,smooth_kernel);
% for i=1:dim
    % b0=ori_normals(:,i);
    % b=H1D*b0;
    % figure;
    % hold on;
    % plot( b, 'r-o' );
    % plot( b0, 'b-x' );
    % hold off;
% end
% pause;
% for i=1:dim
    % % b=normals(:,i);
    % b0=ori_normals(:,i);
    % b=H1D*b0;
    % lambda=0.2;
    % x=tvl2_total_variation_vec_H(b,lambda,A,H1D);
    % figure;
    % hold on;
    % plot( b, 'r-o' );
    % plot( x, 'b-x' );
    % hold off;
% end
% return;

nout= normalOpt(points, normals, lambda, k, false, A , E, H);
nout= normalize_normals(nout);
% plot_1D_cloud(normals,nout);
normals=nout;
% normals = normalOpt(points, normals, lambda, k, true, A , E, H);
% normals = normalize_normals(normals);
draw_points_and_normals(points,normals,'normalOpt');

[points] = optimizePos(points, normals, k, lambda, A, E);
draw_points_and_normals(points,normals,'optimizePos');
