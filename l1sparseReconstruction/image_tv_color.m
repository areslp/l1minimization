function [] = image_tv()
clear all;
close all;

% matlabpool(2)

% I=imread('lena512color.tiff');
I=imread('aierlans.png');
I=im2double(I);
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

[D,E]=image_differencial_matrix(m,n,k);
% D=compute_weight_image(I,E); % reweighted
lambda=1;
xx=tvl1_total_variation_vec_H(b,lambda,D,H);
% xx=total_variation_vec(b,lambda,3,D);
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
