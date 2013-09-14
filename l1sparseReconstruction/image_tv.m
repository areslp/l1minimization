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
figure;
imagesc(H);
pause;
b=H*b;
I=reshape(b,m,n);
% figure;
% imshow(I);
% pause;

[D,E]=image_differencial_matrix(m,n,1);
D=compute_weight_image(I,E); % reweighted
lambda=1;
x=tvl2_total_variation_vec_H(b,lambda,D,H);
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
