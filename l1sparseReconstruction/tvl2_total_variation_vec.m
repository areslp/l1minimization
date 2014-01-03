function [x, history] = group_lasso_feat_split(b, lambda, D)
t_start = tic;
% Global constants and defaults

QUIET    = 0;
MAX_ITER = 200;
RELPRI  = 0;
RELDUAL   = 0;
ABS=1e-5;
REL=1e-5;
% Data preprocessing
[m,n]=size(D);

% number of subsystems
% N = m/ni;
% ADMM solver

rho = 1;

x = zeros(n,1);
z = zeros(m,1); %就是\bar z
u = zeros(m,1);

I=speye(n);
Dt=D';
DtD=Dt*D;
for iter = 1:MAX_ITER
    xold=x;
    % x=T1\(b+rho*D'*(vec(z)+vec(u)));
    x=(I+rho*DtD)\(b+rho*Dt*(z+u));

    Dx=D*x;
    zold = z;
    z=wthresh(Dx-u,'s',lambda/rho);

    u= u + z - Dx;
    % update rho
    % if norm(z-Dx)>0.8*norm(zold-D*xold) % 0<alpha<1, alpha=0.8
        % rho=2*rho; % gamma=2
    % end

    % convergence
    relchg=max(norm(x-xold),norm(z-zold));
    relpri=norm(z-Dx);

    if norm(z-Dx)>0.9*norm(zold-D*xold)
        rho=1.9*rho; 
    end

    fprintf(1,'relchg:%f, relpri:%f\n',relchg,relpri);
    
    if relchg<ABS && relpri<REL
        break;
    end
end

if ~QUIET
    toc(t_start);
end
