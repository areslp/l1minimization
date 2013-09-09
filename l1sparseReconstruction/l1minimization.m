function [t] = l1minimization(A,f,lambda)
% TODO: 这个函数实现有问题，结果很奇怪
% min ||At+f||_1+lambda||t||_2

t_start = tic;
% Global constants and defaults

QUIET    = 0;
MAX_ITER = 100;
RELTOL  = 1e-2;
RELCHG   = 1e-4;
% Data preprocessing

[m,n]=size(A);

% ADMM solver

rho = 1.9;
t = zeros(n,1);
z = zeros(m,1);
u = zeros(m,1);

AtA=A'*A;

for k = 1:MAX_ITER
    % z-update (to be done in parallel)
    % tic
    zold=z;
    z=wthresh(A*t+f-u,'s',1/rho);
    % t=toc;
    % fprintf(1,'update x: %f\n',t);

    % t-update
    % tic
    told = t;
    t=(2*lambda*speye(n)+rho*AtA)\(rho*A'*(z-f+u));
    %toc;
    % fprintf(1,'update z: %f\n',t);

    % tic
    u= u + z - A*t -f;
    % t=toc;
    % fprintf(1,'update u: %f\n',t);

    % rho=min(1e10,1.9*rho);

    % convergence
    % disp('convergence');
    relchg=max(norm(t-told,'fro'),norm(z-zold,'fro'));
    reltol=norm(z-A*t-f);

    fprintf(1,'relchg:%f, reltol:%f\n',relchg,reltol);
    
    if relchg < RELCHG && reltol<RELTOL
        break;
    end

end

if ~QUIET
    toc(t_start);
end
