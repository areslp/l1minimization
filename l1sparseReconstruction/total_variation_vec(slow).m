function [x, history] = group_lasso_feat_split(b, lambda, ni, D, neighbor)
t_start = tic;
% Global constants and defaults

QUIET    = 0;
MAX_ITER = 100;
RELTOL  = 1e-2;
RELCHG   = 1e-4;
% Data preprocessing

n=length(b);
m=n;

% check that ni divides in to n
if (rem(n,ni) ~= 0)
    error('invalid block size');
end
% number of subsystems
N = n/ni
% ADMM solver

rho = 1.9;
alpha = 1.4;    % over-relaxation parameter

x = zeros(n,1);
z = zeros(neighbor*ni,N); %就是\bar z
u = zeros(neighbor*ni,N);
Dix= zeros(neighbor*ni,N);

% fake A
A=eye(neighbor*ni);
[V,Dig] = eig(A'*A);
Dig = diag(Dig);
At = A';


% pre-factor
sum_DitDi=zeros(n);
for i = 1:N
    % (i-1)*neighbor*ni+1
    % i*neighbor*ni
    Di{i}=D((i-1)*neighbor*ni+1:i*neighbor*ni,:); % sub差分矩阵
    Dit{i}=Di{i}';
    DitDi{i}=Dit{i}*Di{i};
    sum_DitDi=sum_DitDi+DitDi{i};
end

for k = 1:MAX_ITER
    % x-update (to be done in parallel)
    xold=x;
    sum_Ditzi=zeros(n,1);
    sum_Ditui=zeros(n,1);
    for i=1:N
        sum_Ditzi=sum_Ditzi+Dit{i}*z(:,i);
        sum_Ditui=sum_Ditui+Dit{i}*u(:,i);
    end
    x=(speye(n)+rho*sum_DitDi)\(b+rho*(sum_Ditzi+sum_Ditui));

    % z-update
    zold = z;
    for i = 1:N
        Dix(:,i)=Di{i}*x;
        xx = x_update(A, Dix(:,i)-u(:,i), lambda/rho, V, Dig);
        z(:,i)=xx;
    end

    % u-update
    for i = 1:N
        u(:,i) = u(:,i) + z(:,i)-Di{i}*x;
    end

    % convergence
    % disp('convergence');
    relchg=max(norm(x-xold,'fro'),norm(z-zold,'fro'));
    tmp=sum(Dix-z);
    reltol=max(tmp(:));

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

q = A'*b;

if (norm(q) <= kappa)
   x = zeros(n,1);
else
    % bisection on t
    lower = 0; upper = 1e10;
    for i = 1:100,
        t = (upper + lower)/2;

        x = V*((V'*q)./(D + t));
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
