function [ Nout ] = cvx_wrapper(Xin, Nin, k, lambda, weighted)
n = size( Xin, 1 );
B=construct_Adj_cvx(Xin,Nin,k,weighted); % B is the adjacent matrix
tn=size(B,1);

% cvx optimization
disp('init cvx solver');
cvx_clear;
% cvx_begin 
    % cvx_solver Mosek 
    % p=ones(tn,1);
    % variable Nout(n,3)
    % variable t(tn,1)
    % minimize(p'*t);
    % subject to
        % sum((Nout-Nin).^2,2)<=lambda^2;
        % sum((B*Nout).^2,2)<=t;
% cvx_end


cvx_begin 
    cvx_solver Mosek 
    variable Nout(n,3);
    T=B*Nout;
    minimize(sum(sum(T.^2,2))+lambda*norm(Nout-Nin,1))
    subject to
        % sum((Nout-Nin).^2,2)<=lambda^2;
cvx_end
