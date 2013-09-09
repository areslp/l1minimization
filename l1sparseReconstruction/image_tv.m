function [] = image_tv()
clear all;
close all;

% I=imread('Lena512.png');
% I=imread('Cameraman.bmp');
I=imread('Barbara.png');
I=im2double(I);
[m,n,k]=size(I);
if k~=1
    I=rgb2gray(I);
end
I=imnoise(I,'gaussian',0,0.01);
imwrite(I,'noise_input.png');
b=reshape(I,m*n,1);
[D]=DiffOper3(m,n,1);
% figure;
% imagesc(D);
% colormap('gray');
% pause;
lambda=0.07;
x=total_variation_vec(b,lambda,1,D);
% x=total_variation(b,lambda,D);
OUT=reshape(x,m,n);
imwrite(OUT,'out.jpg');
figure;
subplot(1,2,1);
imshow(I);
subplot(1,2,2);
imshow(OUT);
disp(['chg:' num2str(norm(I-OUT,'fro'))]);

% real function
function [B Bt BtB] = DiffOper(N)
D = spdiags([-ones(N,1) ones(N,1)], [0 1], N,N+1);
D(:,1) = [];
D(1,1) = 0;
B = [ kron(speye(N),D) ; kron(D,speye(N)) ];
Bt = B';
BtB = Bt*B;

function [B]= DiffOper3(m,n,k)
% except the last column
E0=zeros(m*(n-1)*k,2);
% except the last row
E2=zeros((m-1)*n*k,2);
% the last column
E1=zeros(k*m,2);
% last row
E3=zeros(k*n,2);
index1=1;
index2=1;
for i=1:m
    % each row
    gidx=(i-1)*k*n+1; % 每行第一个像素的起始索引
    assert(gidx>0);
    E0(index1:index1+k*(n-1)-1,1)=linspace(gidx,gidx+k*(n-1)-1,k*(n-1)); % matlab里从i开始num个数这么写：i:i+num-1
    index1=index1+k*(n-1);
    % last column
    E1((i-1)*k+1:(i-1)*k+k,1)=linspace(k*i*n-k+1,k*i*n,k);
    % last row
    if i==m
        E3(:,1)=linspace(gidx,gidx+k*n-1,k*n);
    else
        E2(index2:index2+k*n-1,1)=linspace(gidx,gidx+k*n-1,k*n);
        index2=index2+k*n;
    end
end
E0(:,2)=E0(:,1)+k;
E1(:,2)=E1(:,1)-k;
E2(:,2)=E2(:,1)+k*n;
E3(:,2)=E3(:,1)-k*n;
E=[E0;E1;E2;E3];
rowe=size(E,1);
% assert(rowe==2*m*n*k);
B = sparse( [1:rowe, 1:rowe], [E(:,1); E(:,2)],...
    [ones(rowe, 1), -ones(rowe, 1)], rowe, m*n*k);

% vector valued
% use kdtree to find the nearest neighbor
function [B Bt BtB] = DiffOper2(N)
k=4; % four neighbor
B=sparse([],[],[],k*N*N,N*N,k*N*N);
% generate the point set
points=zeros(N*N,2);
tmp=repmat([1:N]',N,1);
points(:,1)=tmp;
% first column
tmp=[1:N];
tmp=repmat(tmp,N,1);
tmp=tmp(:)';
points(:,2)=tmp;
B=kdtree_adj(points,k);
Bt = B';
BtB = Bt*B;
