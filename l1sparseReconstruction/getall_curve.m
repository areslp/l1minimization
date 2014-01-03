clear all;
close all;

k=6;
lambda=1;

[points,normals]=gen_test();
normals = normalize_normals(normals);
n=size(points,1);
dim=size(points,2);
smooth_kernel=20;
ori_normals=normals;
Asize=n*dim;


% [A,E]=kdtree_adj_vec(points,k,3); % k变化了，需要重新计算
tic
[i,j,v,E]=compute_AE(points,normals,k,dim); % k变化了，需要重新计算
E=E+1;
A=sparse(i+1,j+1,v,n*dim*k,n*dim);
t=toc;
fprintf(1,'compute AE takes:%f\n',t);

% normals=pca_normal(points,k); % TODO: 不是consist的法矢
% draw_points_and_normals(points,normals,'ori');

% N=reshape(normals',dim*n,1); % vectorize
% vc=compute_Dx(A,N,3);
% fprintf(2,'correct Dx: %f\n',vc);

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
draw_points_and_normals(points,normals,'init normal');


tic
d=1;
[i,j,v,E]=compute_AE(points,normals,k,d); % k变化了，需要重新计算
E=E+1;
A=sparse(i+1,j+1,v,n*d*k,n*d);
t=toc;
fprintf(1,'compute AE takes:%f\n',t);
W=compute_weight(points,normals,E);
% D=diag(W)*A;
D=A;
S=D*normals;
obj=sqrt(sum(S.^2,2));
% 截到范围0-1
obj(find(obj>1))=1;
sum(obj)
centers = 0:0.01:1;
counts = hist(obj,centers);
pcts = counts / sum(counts);
figure;
bar(centers,pcts);
ylabel('%');
axis([0 1 0 1]);
axis tight;
set(gcf, 'Color', 'w');

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



% [A,E]=kdtree_adj_vec(points,k,3); % k变化了，需要重新计算
tic
[i,j,v,E]=compute_AE(points,normals,k,dim); % k变化了，需要重新计算
E=E+1;
A=sparse(i+1,j+1,v,n*dim*k,n*dim);
t=toc;
fprintf(1,'compute AE takes:%f\n',t);

% normals=pca_normal(points,k); % TODO: 不是consist的法矢
% draw_points_and_normals(points,normals,'ori');

% N=reshape(normals',dim*n,1); % vectorize
% vc=compute_Dx(A,N,3);
% fprintf(2,'correct Dx: %f\n',vc);

tic
[i,j,v]=compute_WL(points,smooth_kernel);
H=sparse(i+1,j+1,v,Asize,Asize);
t=toc;
fprintf(1,'compute blur kernel takes:%f\n',t);

lambda=10;
nout= normalOpt(points, normals, lambda, k, false, A , E, H);
nout= normalize_normals(nout);
% plot_1D_cloud(ori_normals,nout);
normals=nout;
draw_points_and_normals(points,normals,'normalOpt1');

normals = normalOpt(points, normals, lambda, k, true, A , E, H);
normals = normalize_normals(normals);
draw_points_and_normals(points,normals,'normalOpt2');

tic
d=1;
[i,j,v,E]=compute_AE(points,normals,k,d); % k变化了，需要重新计算
E=E+1;
A=sparse(i+1,j+1,v,n*d*k,n*d);
t=toc;
fprintf(1,'compute AE takes:%f\n',t);
W=compute_weight(points,normals,E);
D=diag(W)*A;
S=D*normals;
obj=sqrt(sum(S.^2,2));
% 截到范围0-1
obj(find(obj>1))=1;
sum(obj)
centers = 0:0.01:1;
counts = hist(obj,centers);
pcts = counts / sum(counts);
figure;
bar(centers,pcts);
ylabel('%');
axis([0 1 0 1]);
axis tight;
set(gcf, 'Color', 'w');


% [A,E]=kdtree_adj_vec(points,k,3); % k变化了，需要重新计算
tic
[i,j,v,E]=compute_AE(points,normals,k,dim); % k变化了，需要重新计算
E=E+1;
A=sparse(i+1,j+1,v,n*dim*k,n*dim);
t=toc;
fprintf(1,'compute AE takes:%f\n',t);

% normals=pca_normal(points,k); % TODO: 不是consist的法矢
% draw_points_and_normals(points,normals,'ori');

% N=reshape(normals',dim*n,1); % vectorize
% vc=compute_Dx(A,N,3);
% fprintf(2,'correct Dx: %f\n',vc);

tic
[i,j,v]=compute_WL(points,smooth_kernel);
H=sparse(i+1,j+1,v,Asize,Asize);
t=toc;
fprintf(1,'compute blur kernel takes:%f\n',t);


lambda=10;
[points] = optimizePos(points, normals, k, lambda, A, E);
draw_points_and_normals(points,normals,'positionOpt');
