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
% L=HtH+rho*DtD;
Htb=Ht*b;
% FH=fftn(full(H));
% FD=fftn(full(D));
% L=norm(FH)*norm(FH)+rho*norm(FD)*norm(FD);
% F=linfactor(L);
for iter = 1:MAX_ITER
    % tic
    xold=x;
    L=HtH+rho*DtD;
    R=Htb+rho*Dt*(z+u);
    % FR=fftn(R);
    % x=ifftn(FR);
    % x=linfactor(F,R);
    x=L\R;
    % toc

    % tic
    Dx=D*x;
    zold = z;
    z=wthresh(Dx-u,'s',lambda/rho);
    % toc

    u= u + z - Dx;

    % update rho
    if norm(z-Dx)>0.8*norm(zold-D*xold) % 0<alpha<1
        rho=2*rho; % gamma=2
    end

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
