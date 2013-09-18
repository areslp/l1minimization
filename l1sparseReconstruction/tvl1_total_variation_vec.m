function [x, history] = group_lasso_feat_split(b, lambda, D)
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
rho1 = 1; % z=Dx
rho2 = 1; % y=x-b

x = zeros(n,1);
z = zeros(m,1); %就是\bar z
y = zeros(n,1);
u1 = zeros(m,1);
u2 = zeros(n,1);

I=speye(n);
Dt=D';
DtD=Dt*D;
L=rho2*I+rho1*DtD;
F=linfactor(L);
for iter = 1:MAX_ITER
    xold=x;
    R=(rho2*(b+y-u2)+rho1*Dt*(z+u1));
    x=linfactor(F,R);
    % x=(rho2*I+rho1*DtD);

    Dx=D*x;
    zold = z;
    z=wthresh(Dx-u1,'s',lambda/rho1); %被拆分成了D_x,D_y,D_t等

    yold = y;
    y=wthresh(x-b+u2,'s',1/(2*rho2));

    u1= u1 + z - Dx;
    u2= u2 + x-b-y;

    relchg=max(norm(x-xold),norm(z-zold));
    relchg=max(relchg,norm(y-yold));
    relpri=max(norm(z-Dx),norm(x-b-y));

    fprintf(1,'relchg:%f, relpri:%f\n',relchg,relpri);
    
    if relchg<ABS && relpri<REL
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
