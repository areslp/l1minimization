function [x, history] = group_lasso_feat_split(b, lambda, ni, D)
t_start = tic;
% Global constants and defaults

QUIET    = 0;
MAX_ITER = 100;
RELPRI  = 0;
RELDUAL   = 0;
ABS=1e-2;
REL=1e-4;
% Data preprocessing
[m,n]=size(D);

% number of subsystems
N = m/ni;
% ADMM solver

rho = 10;
mu=10;
alpha = 1.4;    % over-relaxation parameter

x = zeros(n,1);
z = zeros(m,1); %就是\bar z
u = zeros(m,1);

bf=norm(b,'fro');
% T1=speye(n)+rho*D'*D;

for k = 1:MAX_ITER
    xold=x;
    % x=T1\(b+rho*D'*(vec(z)+vec(u)));
    x=(speye(n)+rho*D'*D)\(b+rho*D'*(vec(z)+vec(u)));

    Dx=D*x;
    zold = z;
    z=wthresh(Dx-u,'s',lambda/rho);

    r=z-Dx;
    s=rho*(-D)'*(z-zold);

    u= u + z - Dx;

    if norm(r)>mu*norm(s)
        rho=2*rho;
    elseif norm(s)>mu*norm(r)
        rho=rho/2;
    else
        % do nothing
    end

    % convergence
    RELPRI=sqrt(m)*ABS+REL*max(norm(-D*x,2),max(norm(z),0));
    RELDUAL=sqrt(n)*ABS+REL*norm((-D)'*u*rho); % u=y/rho

    fprintf(1,'relpri:%f, reldual:%f, norm(r):%f, norm(s):%f\n',RELPRI,RELDUAL,norm(r),norm(s));
    
    if norm(r)< RELPRI && norm(s) < RELDUAL
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
