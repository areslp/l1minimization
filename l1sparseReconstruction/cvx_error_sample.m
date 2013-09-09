load vars.mat;
cvx_clear;
% cvx optimization
cvx_begin 
    % cvx_solver sedumi
    cvx_solver gurobi 
    p=ones(kn,1);
    variable Nout(n,3)
    variable t(kn,1)
    minimize(p'*t);
    subject to
        sum((Nout-Nin).^2,2)<=lambda;
        sum((B*Nout).^2,2)<=t;
        % for i=1:kn
            % Ai=B(i,:);
            % norm(Ai*Nout)<=t(i);
        % end
cvx_end
