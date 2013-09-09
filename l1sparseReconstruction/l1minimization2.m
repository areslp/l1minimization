function [x]=l1minimization2(W,D,b,lambda)
% min_x 1/2||x-b||_2^2+\lambda \sum||w_i'*D_i*x||_1

% W kx1 cell, each cell 3x1, kx3
% D 3kxn

t_start = tic;
% Global constants and defaults

QUIET    = 0;
MAX_ITER = 100;
RELTOL  = 1e-2;
RELCHG   = 1e-4;

% Data preprocessing
k=size(W,1);
n=size(D,2);
% W kx3k
% w1 0  0 ...
% 0  w2 0 ...
% 0  0  w3 ...
% .
% .
% .
% 生成[1 1 1 0 0....]这样的行向量
% c=mat2cell(W,ones(k,1),3);
% disp('mat2cell end');
% W=blkdiag(c{:});
% disp('blkdiag end');

W=W'; % 3xk

% 这个方法绝妙
[r,c] = size(W);
i     = 1:numel(W);
j     = repmat(1:c,r,1);
W     = sparse(i',j(:),W(:)); % j(:)是按列优先排成的一个向量


disp(['size W:' num2str(size(W))]);
disp(['size D:' num2str(size(D))]);

W=W';

A=W*D;
size(A)

% A的计算方法
% tic
% for i=1:k
    % wi=W(i,:); % 1x3
    % Di=D(3*(i-1)+1:3*i,:); % 3xn
    % Ai=wi*Di;
    % A(i,:)=Ai;
% end
% toc

AtA=A'*A;

% ADMM solver

rho = 1.9;
x = zeros(n,1);
z = zeros(k,1);
u = zeros(k,1);

bf=norm(b,'fro');
% T1=speye(n)+rho*AtA;

for k = 1:MAX_ITER
    zold=z;
    z=wthresh(A*x-u,'s',lambda/rho);

    xold = x;
    % x=T1\(b+(rho*A'*(z+u)));
    x=(speye(n)+rho*AtA)\(b+(rho*A'*(z+u)));

    u= u + z - A*x;
    % rho=min(1e10,1.9*rho);

    % convergence
    % 第一次迭代就全部都=0，这是怎么回事？
    relchg=max(norm(x-xold,'fro'),norm(z-zold,'fro'))/bf;
    reltol=norm(z-A*x,'fro')/bf;

    fprintf(1,'relchg:%f, reltol:%f\n',relchg,reltol);
    
    if relchg < RELCHG && reltol<RELTOL
        break;
    end
end

if ~QUIET
    toc(t_start);
end
