function [points,normals] = gen_test()
step=0.01;
nl = 20;

x2 = 0:step:1;
x1 = -1:step:0;
y2 = -x1;
y1 = x2;
% x3=zeros(1,length(x1));
% y3=1:step:2;
y1 = awgn(y1,nl,'measured');
y2 = awgn(y2,nl,'measured');

p1=[x1;y1;zeros(1,length(x1))]';
p2=[x2;y2;zeros(1,length(x2))]';
% p3=[x3;y3;zeros(1,length(x3))]';
% p=[p1;p2;p3];
p=[p1;p2];
[row,col] = size(p);
normals=zeros(row,3);
% normals(1:row/3,:)=repmat([-1 1 0],row/3,1);
% normals(row/3+1:2*row/3,:)=repmat([1 1 0],row/3,1);
% normals(2*row/3+1:end,:)=repmat([1 0 0],row/3,1);

normals(1:row/2,:)=repmat([-1 1 0],row/2,1);
normals(row/2+1:end,:)=repmat([1 1 0],row/2,1);

points=p;
[p,ia,ic]=unique(points,'rows');
% p=points(ia,:);
% points=p(ic,:);
normals=normals(ia,:);
points=p;
