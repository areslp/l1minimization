function [] = image_tv()
clear all;
close all;

% matlabpool(2)

I=imread('flowers.png');
I=im2double(I);
I0=I;
b0=color2vector(I);
[m,n,k]=size(I);
b=color2vector(I);

[D,E]=image_differencial_matrix(m,n,k);

iter=2;
DETAILS=cell(iter,1);
lambdas=[0.05 0.1 0.2 0.4 0.8];

for i=1:iter
    b=color2vector(I);
    lambda=lambdas(i);
    xx=total_variation_vec(b,lambda,3,D);
    % xx=(xx-min(xx))/(max(xx)-min(xx));
    OUT=vec2color(xx,m,n);
    BASE=OUT;
    DETAILS{i}=I-OUT;
    I=OUT;
end

weights=[2 2 1 1 1];
ENHANCE=BASE;
for i=1:iter
    ENHANCE=ENHANCE+weights(i)*DETAILS{i};
end
% ENHANCE=(ENHANCE-min(ENHANCE(:)))/(max(ENHANCE(:))-min(ENHANCE(:)));

figure;
subplot(1,2,1);
imshow(I0);
subplot(1,2,2);
imshow(ENHANCE);
norm(color2vector(I0)-color2vector(ENHANCE))

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
