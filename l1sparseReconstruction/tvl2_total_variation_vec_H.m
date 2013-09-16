function [x, history] = group_lasso_feat_split(b, lambda, D, H)
t_start = tic;
% Global constants and defaults

QUIET    = 0;
MAX_ITER = 5000;
RELPRI  = 0;
RELDUAL   = 0;
ABS=1e-5;
REL=1e-5;
% Data preprocessing
[m,n]=size(D);

rho = 1;

x = zeros(n,1);
z = zeros(m,1); %就是\bar z
u = zeros(m,1);

Dt=D';
DtD=Dt*D;
Ht=H';
HtH=Ht*H;
L=HtH+rho*DtD;
Htb=Ht*b;
F=linfactor(L);
for iter = 1:MAX_ITER
    % tic
    xold=x;
    R=Htb+rho*Dt*(z+u);
    x=linfactor(F,R);
    % toc

    % tic
    Dx=D*x;
    zold = z;
    z=wthresh(Dx-u,'s',lambda/rho);
    % toc

    u= u + z - D*x;

    % convergence
    relchg=max(norm(x-xold),norm(z-zold));
    relpri=norm(z-Dx);

    if mod(iter,100)==0
        fprintf(1,'iter:%d, relchg:%f, relpri:%f\n',iter,relchg,relpri);
    end
    
    if relchg<ABS && relpri<REL
        break;
    end
end

if ~QUIET
    toc(t_start);
end
