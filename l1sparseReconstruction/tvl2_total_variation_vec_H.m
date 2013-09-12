function [x, history] = group_lasso_feat_split(b, lambda, D, H)
t_start = tic;
% Global constants and defaults

QUIET    = 0;
MAX_ITER = 1000;
RELPRI  = 0;
RELDUAL   = 0;
ABS=1e-3;
REL=1e-3;
% Data preprocessing
[m,n]=size(D);

rho = 1;

x = zeros(n,1);
z = zeros(m,1); %就是\bar z
u = zeros(m,1);

DtD=D'*D;
for iter = 1:MAX_ITER
    xold=x;
    x=(H'*H+rho*DtD)\(H'*b+rho*D'*(z+u));

    Dx=D*x;
    zold = z;
    z=wthresh(D*x-u,'s',lambda/rho);

    u= u + z - D*x;

    % convergence
    relchg=max(norm(x-xold),norm(z-zold));
    relpri=norm(z-Dx);

    fprintf(1,'relchg:%f, relpri:%f\n',relchg,relpri);
    
    if relchg<ABS && relpri<REL
        break;
    end
end

if ~QUIET
    toc(t_start);
end
