function [B,E]= DiffOper3(m,n,k)
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
