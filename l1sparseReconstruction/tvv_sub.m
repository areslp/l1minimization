function [] = tvv_sub(Ai,v,Di,lambda)
% min 1/2||Aixi-v||_2^2+lambda||Dixi||_2
QUIET    = 0;
MAX_ITER = 100;
RELTOL  = 1e-2;
ABSTOL   = 1e-4;

rho = RHO;
alpha = ALPHA;    % over-relaxation parameter

[~,ni]=size(Ai);

x = zeros(ni,N);
z = zeros(m,1); %就是\bar z
u = zeros(m,1);


        xx = x_update(Ai, Aixi(:,i) + z - Axbar - u, lambda/rho, V{i}, D{i}, E{i});


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

end
