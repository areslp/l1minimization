function [x, history] = group_lasso_feat_split(b, lambda, dim, D)
t_start = tic;
% Global constants and defaults

QUIET    = 0;
MAX_ITER = 200;
RELTOL  = 1e-5;
RELCHG   = 1e-5;
% Data preprocessing
[m,n]=size(D);
% number of subsystems
N = m/dim;
% ADMM solver
rho = 20;

bf=norm(b,'fro');

x = zeros(n,1);
z = zeros(m,1);
u = zeros(m,1);

I=speye(n);
Dt=D';
DtD=D'*D;
L=I+rho*DtD;
F=linfactor(L);
for iter = 1:MAX_ITER
    % x-update (to be done in parallel)
    xold=x;
    % L=I+rho*DtD;
    R=b+rho*Dt*(z+u);
    % x=L\R;
    x=linfactor(F,R);

    % z-update
    Dx=D*x;
    zold = z;
    Dxu=Dx-u;
    z=z_loop(Dxu,lambda/rho,dim);

    u= u + z - Dx;

    % convergence
    relchg=max(norm(x-xold),norm(z-zold))/bf;
    reltol=norm(Dx-z)/bf;

    fprintf(1,'iter:%d, relchg:%f, reltol:%f\n',iter,relchg,reltol);
    
    if relchg < RELCHG && reltol<RELTOL
        break;
    end

    % rho=1.3*rho;
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

function [Xv] = vec(X)
% This operator stacks the columns of the matrix X into a 
% single column vector Xv.
[a,b] = size(X);
Xv = reshape(X,a*b,1);
