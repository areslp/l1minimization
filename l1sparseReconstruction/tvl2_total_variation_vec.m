function [x, history] = group_lasso_feat_split(b, lambda, D)
t_start = tic;
% Global constants and defaults

QUIET    = 0;
MAX_ITER = 1000;
RELPRI  = 0;
RELDUAL   = 0;
ABS=1e-2;
REL=1e-2;
% Data preprocessing
[m,n]=size(D);

% number of subsystems
% N = m/ni;
% ADMM solver

rho = 10;
mu=10;

x = zeros(n,1);
z = zeros(m,1); %就是\bar z
u = zeros(m,1);

% T1=speye(n)+rho*D'*D;
I=speye(n);
DtD=D'*D;
for iter = 1:MAX_ITER
    xold=x;
    % x=T1\(b+rho*D'*(vec(z)+vec(u)));
    x=(I+rho*DtD)\(b+rho*D'*(z+u));

    Dx=D*x;
    zold = z;
    z=wthresh(Dx-u,'s',lambda/rho);

    u= u + z - Dx;

    % convergence
    relchg=max(norm(x-xold),norm(z-zold));
    relpri=norm(z-Dx);

    % if norm(z-Dx)>0.9*norm(zold-D*xold)
        % rho=1.9*rho; 
    % end

    fprintf(1,'relchg:%f, relpri:%f\n',relchg,relpri);
    
    if relchg<ABS && relpri<REL
        break;
    end
end

if ~QUIET
    toc(t_start);
end
