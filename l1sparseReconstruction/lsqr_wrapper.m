function [NO] = lsqr_wrapper(B,x)

disp('lsqr');
atol   = 1.0e-6;
btol   = 1.0e-6;
conlim = 1.0e+10;
itnlim = 10*size(B,2);
show   = 0;
damp   = 0;

[ NO, istop, itn, r1norm, r2norm, Anorm, Acond, Arnorm, xnorm, var ] ...
  = lsqrSOL( size(B,1), size(B,2), B, x, damp, atol, btol, conlim, itnlim, show );
