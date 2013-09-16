function [x] = constraint_least_square(D,y,z,lambda)
% min 1/2||x-y||_2^2 s.t. z=Dx
RELCHG=1e-5;
RELTOL=1e-5;
rho=1;
[m,n]=size(D);
x=zeros(n,1);
u=zeros(m,1);
MAX_ITER=1000;
Dt=D';
DtD=Dt*D;
I=speye(n,n);
L=lambda*I+rho*DtD;
F=linfactor(L);
R=lambda*y+rho*Dt*(z+u);
x=linfactor(F,R);
return;

for iter=1:MAX_ITER
    % update x
    xold=x;
    R=y+rho*Dt*(z+u);
    x=linfactor(F,R);

    Dx=D*x;
    % update u
    u=u+z-Dx;

    % convergence
    relchg=max(norm(x-xold));
    reltol=norm(Dx-z);

    fprintf(1,'iter:%d, relchg:%f, reltol:%f\n',iter,relchg,reltol);
    
    if relchg < RELCHG && reltol<RELTOL
        break;
    end
end
