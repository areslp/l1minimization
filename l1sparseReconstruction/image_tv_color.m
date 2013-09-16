function [] = image_tv()
clear all;
close all;

% matlabpool(2)

% I=imread('lena512color.tiff');
I=imread('aierlans.png');
I=im2double(I);
b0=color2vector(I);
[m,n,k]=size(I);

sigma=3;
% H       = fspecial('gaussian', [9 9], 2);
% G       = imfilter(I, H, 'circular');
% figure;
% imshow(G);
% add noise
I       = imnoise(I, 'gaussian', 0.05);

% image to vector
b=color2vector(I);

% 分通道进行模糊
% H=generate_PSF_matrix(sigma,m*n);
% for i=1:k
    % I(:,:,i)=reshape(H*reshape(I(:,:,i),m*n,1),m,n); % row blur
    % I(:,:,i)=reshape(H*reshape(I(:,:,i)',n*m,1),n,m)'; % col blur
% end
% imshow(I);

% 直接进行模糊
H=generate_PSF_matrix_vec(m,n,k);
b=H*b;
I=vec2color(b,m,n);
% figure;
% imshow(I);
% pause;

% [D,E]=image_differencial_matrix(m,n,k);
% D=compute_weight_image(I,E); % reweighted

% 把图像看成点云，测试compute_AE是否正确
points=fake_points_from_image(I);
pn=size(points,1);
dim=3;
k=12;
tic
[i,j,v,E]=compute_AE(points,points,k,dim); % k变化了，需要重新计算
E=E+1;
A=sparse(i+1,j+1,v,pn*dim*k,pn*dim);
AA = sparse( [1:1:dim*k*pn, 1:1:dim*k*pn], [E(:,1); E(:,2)],...
    [ones(dim*k*pn, 1), -ones(dim*k*pn, 1)], dim*k*pn, dim*pn);
t=toc;
fprintf(1,'compute AE takes:%f\n',t);
% [A,E]=kdtree_adj_vec(points,k,dim);
% size(A)
% size(H)
D=AA;

vc=compute_Dx(D,b0,3);
v0=compute_Dx(D,b,3);
lambda=0.2;
% xx=tvl2_total_variation_vec_H(b,lambda,D,H);
xx=total_variation_vec_H(b,lambda,3,D,H);
% xx=total_variation_vec(b,lambda,3,D);
v=compute_Dx(D,xx,3);
fprintf(2,'correct Dx is %f\n',vc);
fprintf(2,'init Dx is %f\n',v0);
fprintf(2,'L1 Dx is %f\n',v);
xx=(xx-min(xx))/(max(xx)-min(xx));
OUT=vec2color(xx,m,n);
plot_1D(I,OUT,floor(m/2));
figure;
subplot(1,2,1);
imshow(I);
subplot(1,2,2);
imshow(OUT);

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
