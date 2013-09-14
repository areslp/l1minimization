function [x, history] = group_lasso_feat_split(b, lambda, D, H)
% min 1/2|x-b|_1+lambda|Dx|_1
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

% ADMM solver
rho1 = 10; % z=Dx
rho2 = 10; % y=x-b
mu=10;
alpha = 1.6;    % over-relaxation parameter

x = zeros(n,1);
z = zeros(m,1); %就是\bar z
y = zeros(n,1);
u1 = zeros(m,1);
u2 = zeros(n,1);

I=speye(n);
DtD=D'*D;
HtH=H'*H;
L=(rho2*HtH+rho1*DtD);
F=linfactor(L);
for iter = 1:MAX_ITER
    xold=x;
    R=(rho2*H'*(b+y-u2)+rho1*D'*(z+u1));
    x=linfactor(F,R);

    Dx=D*x;
    zold = z;
    z=wthresh(Dx-u1,'s',lambda/rho1); %被拆分成了D_x,D_y,D_t等

    yold = y;
    y=wthresh(H*x-b+u2,'s',1/(2*rho2));

    u1= u1 + z - Dx;
    u2= u2 + H*x-b-y;

    relchg=max(norm(x-xold),norm(z-zold));
    relchg=max(relchg,norm(y-yold));
    relpri=max(norm(z-Dx),norm(H*x-b-y));

    fprintf(1,'relchg:%f, relpri:%f\n',relchg,relpri);
    
    if relchg<ABS && relpri<REL
        break;
    end

    % acc convergence
    % if norm(z-Dx)/norm(zold-D*xold)>0.7
        % rho1=2*rho1;
    % end
    % if norm(H*x-b-y)/norm(xold-b-yold)>0.7
        % rho2=2*rho2;
    % end 
end

if ~QUIET
    toc(t_start);
end
