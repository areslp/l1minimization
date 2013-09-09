clear;
close all;

data=load('data/double-torus2.xyzn');
% data=load('data/double-torus1.xyzn');
P=data(:,1:3);
N=data(:,4:6);

% parameters
k=6;

% normalize N
N = N./repmat(sqrt(sum(N.^2,2)),1,size(N,2));
n=size(P,1);
tic
Nnew=reconstruct_normal(P,N,k,false);
toc
ff=fopen('data/out.xyzn','w');
for i=1:n
    fprintf(ff,'%f %f %f %f %f %f\n',P(i,1),P(i,2),P(i,3),Nnew(i,1),Nnew(i,2),Nnew(i,3));
end
fclose(ff);
