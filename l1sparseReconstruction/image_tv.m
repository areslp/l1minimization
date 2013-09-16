function [] = image_tv()
clear all;
close all;

% I=imread('Lena512.png');
I=imread('aierlans.png');
% I=imread('Cameraman.bmp');
% I=imread('Barbara.png');
I=im2double(I);
[m,n,k]=size(I);
if k~=1
    I=rgb2gray(I);
    k=1;
end

% H=fspecial('gaussian', [9 9], 2);
% G=imfilter(I, H, 'circular');
% figure;
% imshow(G);

% I = I + 0.05*randn(size(I));
b=reshape(I,m*n,1);

H=generate_PSF_matrix2(m,n);
% figure;
% imagesc(H);
% pause;
b=H*b;
I=reshape(b,m,n);
% figure;
% imshow(I);
% pause;

% 把图像看成点云，测试compute_AE是否正确
points=fake_points_from_image(I);
pn=size(points,1);
dim=1;
k=8;
tic
[i,j,v,E]=compute_AE(points,points,k,1); % k变化了，需要重新计算
E=E+1;
A=sparse(i+1,j+1,v,pn*dim*k,pn*dim);
t=toc;
fprintf(1,'compute AE takes:%f\n',t);
% [A,E]=kdtree_adj_vec(points,k,1);
D=A;
% figure;
% imagesc(D);
% pause;


% [D,E]=image_differencial_matrix(m,n,1);
% D=compute_weight_image(I,E); % reweighted
lambda=0.2;
v=compute_Dx(D,b,1);
fprintf(1,'init Dx is %f\n',v);
% x=tvl2_total_variation_vec_H(b,lambda,D,H);
x=total_variation_vec_H(b,lambda,1,D,H);
v=compute_Dx(D,x,1);
fprintf(1,'L1 Dx is %f\n',v);
x=(x-min(x))/(max(x)-min(x));
OUT=reshape(x,m,n);
imwrite(OUT,'out.jpg');
figure;
subplot(1,2,1);
imshow(I);
subplot(1,2,2);
imshow(OUT);
disp(['chg:' num2str(norm(I-OUT,'fro'))]);
plot_1D(I,OUT,floor(m/2));
