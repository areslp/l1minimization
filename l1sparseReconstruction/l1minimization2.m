function [x]=l1minimization2(W,D,b,lambda)
% min_x 1/2||x-b||_2^2+\lambda \sum||w_i'*D_i*x||_1

% W knx3
% D 3kxn

t_start = tic;
% Global constants and defaults

QUIET    = 0;
MAX_ITER = 100;
RELTOL  = 1e-2;
RELCHG   = 1e-4;

% Data preprocessing
kn=size(W,1);
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

W=W'; % 3xkn

% 这个方法绝妙
[r,c] = size(W);
i     = 1:numel(W);
j     = repmat(1:c,r,1);
W     = sparse(i',j(:),W(:)); % j(:)是按列优先排成的一个向量


disp(['size W:' num2str(size(W))]);
disp(['size D:' num2str(size(D))]);

W=W';

A=W*D;
disp(['size A:' num2str(size(A))]);
At=A';
% A的计算方法
% tic
% for i=1:k
    % wi=W(i,:); % 1x3
    % Di=D(3*(i-1)+1:3*i,:); % 3xn
    % Ai=wi*Di;
    % A(i,:)=Ai;
% end
% toc

AtA=At*A;

% ADMM solver

rho = 1;
x = zeros(n,1);
z = zeros(kn,1);
u = zeros(kn,1);

bf=norm(b,'fro');
L=speye(n)+rho*AtA;
F=linfactor(L);
for iter = 1:MAX_ITER
    zold=z;
    z=wthresh(A*x-u,'s',lambda/rho);

    xold = x;
    % x=T1\(b+(rho*A'*(z+u)));
    % x=(speye(n)+rho*AtA)\(b+(rho*At*(z+u)));
    R=b+rho*At*(z+u);
    % x=L\R;
    x=linfactor(F,R);

    u= u + z - A*x;

    % convergence
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
