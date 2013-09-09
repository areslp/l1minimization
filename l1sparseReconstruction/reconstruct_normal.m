function [Nout] = reconstruct_normal(P, Nin, k, weighted)
% P n*3
% Nin n*3
B=construct_Adj(P,Nin,k,true);
n=size(P,1);
fprintf(1,'point size is %d\n',n);

N=B*Nin;

% cvx optimization
tau=0.1;
cvx_begin 
    variable X(3*k*n,1)
    Obj=sum_square(X,2);
    minimize(sum(Obj))
    subject to
        % for i=1:n
            % norm(Nout(i,:)-Nin(i,:))<=tau;
        % end
        sum((X - Q).^2,2) <= tau;
cvx_end


% disp('lsqr');
% atol   = 1.0e-6;
% btol   = 1.0e-6;
% conlim = 1.0e+10;
% itnlim = 10*size(A,2);
% show   = 0;
% damp   = 0;

% [ Nout, istop, itn, r1norm, r2norm, Anorm, Acond, Arnorm, xnorm, var ] ...
    % = lsqrSOL( size(A,1), size(A,2), A, X, damp, atol, btol, conlim, itnlim, show );

Nout=A\X;
