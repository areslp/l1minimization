function [Xout] = optimizePos(points, normals, k, lambda, A ,E)
n = size(points,1);
dim=size(normals,2);
%% 计算N_avg,E中两点向量均值
Navg = zeros(k * n, dim);
for i = 1:k*n
    idx1=E((i-1)*dim+1,1);
    idx2=E((i-1)*dim+1,2);
    if rem(idx1,dim)==0
        idx1=floor(idx1/dim)-1;
    else
        idx1=floor(idx1/dim);
    end
    idx1=idx1+1;
    if rem(idx2,dim)==0
        idx2=floor(idx2/dim)-1;
    else
        idx2=floor(idx2/dim);
    end
    idx2=idx2+1;
    n1 = normals(idx1,:);
    n2 = normals(idx2,:);
    Navg(i, :) = (n1 + n2) ./ 2; % TODO: 这样加会有些问题，比如正好相反的两个法矢相加为0
end
fprintf('finish calculating Navg %dx%d\n', size(Navg,1), size(Navg,2));

D=A;
b=reshape(points',dim*n,1); % vectorize
x=l1minimization2(Navg,D,b,lambda); % lambda越大，变化越小

%% 得出Xout
Xout=reshape(x,dim,n)';
end
