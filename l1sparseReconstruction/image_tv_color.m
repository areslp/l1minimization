function [] = image_tv()
clear all;
close all;

% matlabpool(2)

% I=imread('color_small.png');
I=imread('color.jpg');
% I=imread('1s.png');
% I=imread('1.jpg');
I=im2double(I);

I_noisy = double(I) + 0.1.*randn(size(I));
% I_noisy = imnoise(I,'gaussian', 0, 0.1);

[m,n,k]=size(I);
assert(m==n);
b=color2vector(I_noisy);
[D]=DiffOper3(m,n);
size(D)
% figure;
% imagesc(D);
% colormap('gray');
% pause;
lambda=0.08;

% METHOD1
% x=group_lasso_wrapper(D, speye(m*n*k*4), D*b, lambda);
% x=(x-min(x))/(max(x)-min(x));
% [xx]=lsqr_wrapper(D,x);
% xx=(xx-min(xx))/(max(xx)-min(xx));
% assert(length(xx)==m*n*k);
% OUT=vec2color(xx,m,n);

% METHOD2
xx=total_variation_vec(b,lambda,3,D,2);
xx=(xx-min(xx))/(max(xx)-min(xx));
OUT=vec2color(xx,m,n);

norm(xx-b)

imwrite(OUT,'out.jpg');
figure;
subplot(1,3,1);
imshow(I);
subplot(1,3,2);
imshow(I_noisy);
subplot(1,3,3);
imshow(OUT);

% matlabpool close;

function [B]= DiffOper3(m,n)
% except the last column
E0=zeros(m*(n-1)*3,2);
% except the last row
E2=zeros((m-1)*n*3,2);
% the last column
E1=zeros(3*m,2);
% last row
E3=zeros(3*n,2);
index1=1;
index2=1;
for i=1:m
    % each row
    gidx=(i-1)*3*n+1;
    assert(gidx>0);
    E0(index1:index1+3*(n-1)-1,1)=linspace(gidx,gidx+3*(n-1)-1,3*(n-1));
    index1=index1+3*(n-1);

    % last column
    E1((i-1)*3+1:(i-1)*3+3,1)=linspace(3*i*n-2,3*i*n,3);

    % last row
    if i==m
        E3(:,1)=linspace(gidx,gidx+3*n-1,3*n);
    else
        E2(index2:index2+3*n-1,1)=linspace(gidx,gidx+3*n-1,3*n);
        index2=index2+3*n;
    end
end
E0(:,2)=E0(:,1)+3;
E1(:,2)=E1(:,1)-3;
E2(:,2)=E2(:,1)+3*n;
E3(:,2)=E3(:,1)-3*n;
E=[E0;E1;E2;E3];
rowe=size(E,1);
assert(rowe==2*m*n*3);
B = sparse( [1:rowe, 1:rowe], [E(:,1); E(:,2)],...
    [ones(rowe, 1), -ones(rowe, 1)], rowe, m*n*3);

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

