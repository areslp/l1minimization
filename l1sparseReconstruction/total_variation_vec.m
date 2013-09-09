function [x, history] = group_lasso_feat_split(b, lambda, ni, D)
t_start = tic;
% Global constants and defaults

QUIET    = 0;
MAX_ITER = 100;
RELTOL  = 1e-2;
RELCHG   = 1e-4;
% Data preprocessing
[m,n]=size(D);

% check that ni divides in to n
if (rem(n,ni) ~= 0)
    error('invalid block size');
end
% number of subsystems
N = m/ni;
% ADMM solver

rho = 1.9;
alpha = 1.4;    % over-relaxation parameter

x = zeros(n,1);
z = zeros(ni,N); %就是\bar z
u = zeros(ni,N);

bf=norm(b,'fro');
% T1=speye(n)+rho*D'*D;

for k = 1:MAX_ITER
    % x-update (to be done in parallel)
    % tic
    xold=x;
    % x=T1\(b+rho*D'*(vec(z)+vec(u)));
    x=(speye(n)+rho*D'*D)\(b+rho*D'*(vec(z)+vec(u)));
    % t=toc;
    % fprintf(1,'update x: %f\n',t);

    % z-update
    Dx=D*x;
    Dx=reshape(Dx,ni,N);
    % tic
    zold = z;
    % 显示解
    % for i = 1:N
        % z(:,i)=vst(lambda/rho, Dx(:,i)-u(:,i));
    % end
    z=z_loop(Dx-u,lambda/rho);
    % norm(z) 
    % t=toc;
    % fprintf(1,'update z: %f\n',t);

    % tic
    u= u + z - Dx;
    % t=toc;
    % fprintf(1,'update u: %f\n',t);

    % rho=min(1e10,1.9*rho);

    % convergence
    % disp('convergence');
    relchg=max(norm(x-xold,'fro'),norm(z-zold,'fro'))/bf;
    reltol=norm(Dx-z,'fro')/bf;

    fprintf(1,'relchg:%f, reltol:%f\n',relchg,reltol);
    
    if relchg < RELCHG && reltol<RELTOL
        break;
    end

end

if ~QUIET
    toc(t_start);
end


% P69 subproblem
function x = x_update(A, b, kappa, V, D)
[m,n] = size(A);

% tic
q = A'*b;
% toc

if (norm(q) <= kappa)
   x = zeros(n,1);
else
    % bisection on t
    lower = 0; upper = 1e10;
    for i = 1:100,
        t = (upper + lower)/2;

        x = V*((V'*q)./(D + t));
        % size(V) nxn
        % size(q) nx1
        % size(D) nx1
        % size(x) nx1
        % pause;
        if t > kappa/norm(x),
            upper = t;
        else
            lower = t;
        end
        if (upper - lower <= 1e-6)
            break;
        end
    end
end

% vector soft thresholding operator 
function [a] = vst(k,v)
    if norm(v)==0
        a=0;
        return;
    end
    t = (1-k/norm(v));
    t = max(t,0);
    a = t*v;
