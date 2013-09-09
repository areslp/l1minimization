function [x] = total_variation(b, lambda, D)
% total_variation  Solve total variation minimization via ADMM
%
% [x, history] = total_variation(b, lambda, rho, alpha)
% 
% Solves the following problem via ADMM:
% 
%   minimize  (1/2)||x - b||_2^2 + lambda * ||Dx||_1
%
% where b in R^n.
%
% The solution is returned in the vector x.

t_start = tic;

%% Global constants and defaults

MAX_ITER = 1000;
ABSTOL   = 1e-2;
RELTOL   = 1e-4;
mu = 1.9;
max_mu = 1e10;
rho=1.9;

%% Data preprocessing

n = length(b);
[k,n]=size(D);

%% ADMM solver

x = zeros(n,1);
z = zeros(k,1);
y = zeros(k,1);

I = speye(n);
DtD = D'*D;

bf=norm(b);

for k = 1:MAX_ITER

    % x-update
    xold=x;
    x = (I + mu*DtD) \ (b + D'*(mu*z+y));

    % z-update with relaxation
    zold = z;
    Dx=D*x;
    z=wthresh(Dx-y/mu,'s',lambda/mu); 

    % y-update
    y = y + mu*(z-Dx);
    % mu = min(max_mu,rho*mu);

    % check convergence
    reltol = max(norm(x-xold,'fro')/bf, norm(z-zold,'fro')/bf);
    abstol = norm(z-Dx,'fro')/bf;

    fprintf(1,'reltol: %f, abstol: %f\n',reltol,abstol);
    
    if reltol < RELTOL && abstol < ABSTOL
        break; 
    end

end

toc(t_start);

end

function obj = objective(b, lambda, D, x, z)
    obj = .5*norm(x - b)^2 + lambda*norm(z,1);
end
