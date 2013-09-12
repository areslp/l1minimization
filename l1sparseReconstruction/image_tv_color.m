function [] = image_tv()
clear all;
close all;

% matlabpool(2)

I=imread('lena512color.tiff');
% I=imread('aierlan.png');
I=im2double(I);

[rows cols frames] = size(I);
% H       = fspecial('gaussian', [9 9], 2);
% I       = imfilter(I, H, 'circular');
% I       = imnoise(I, 'salt & pepper', 0.05);
I       = imnoise(I, 'gaussian', 0.05);

[m,n,k]=size(I);
b=color2vector(I);
[D,E]=image_differencial_matrix(m,n,3);
% D=compute_weight_image(I,E); % reweighted
lambda=0.2;
% METHOD1
% x=group_lasso_wrapper(D, speye(m*n*k*4), D*b, lambda);
% x=(x-min(x))/(max(x)-min(x));
% [xx]=lsqr_wrapper(D,x);
% xx=(xx-min(xx))/(max(xx)-min(xx));
% assert(length(xx)==m*n*k);
% OUT=vec2color(xx,m,n);

% METHOD2
% xx=tvl2_total_variation_vec(b,lambda,D);
% xx=tvl1_total_variation_vec(b,lambda,D);
xx=total_variation_vec(b,lambda,3,D);
xx=(xx-min(xx))/(max(xx)-min(xx));
OUT=vec2color(xx,m,n);

plot_1D(I,OUT,floor(m/2));
norm(xx-b)

imwrite(OUT,'out.jpg');
figure;
subplot(1,2,1);
imshow(I);
subplot(1,2,2);
imshow(OUT);

% vector valued
% use kdtree to find the nearest neighbor
function [B] = DiffOper2(N)
k=4; % four neighbor
% B=sparse([],[],[],k*N*N,N*N,k*N*N);
% generate the point set
points=zeros(N*N,2);
tmp=repmat([1:N]',N,1);
points(:,1)=tmp;
% first column
tmp=[1:N];
tmp=repmat(tmp,N,1);
tmp=tmp(:)';
points(:,2)=tmp;
B=kdtree_adj_vec(points,k,3);

function [vec] = color2vector(I)
[m,n,k]=size(I);
vec=zeros(m*n*k,1);
% 先行后列
idx=1;
for i=1:m
    for j=1:n
        for l=1:k
            vec(idx)=I(i,j,l); 
            idx=idx+1;
        end
    end
end

function [I] = vec2color(vec,m,n)
I=zeros(m,n,3);
% 先行后列
idx=1;
for i=1:m
    for j=1:n
        for k=1:3
            I(i,j,k)=vec(idx); 
            idx=idx+1;
        end
    end
end
