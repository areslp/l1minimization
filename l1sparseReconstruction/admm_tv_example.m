%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%% Sample signal denoising & deblurring using ADMM                 %%%
%%% code written by James Gregson (james.gregson@gmail.com), 2013       %%%
%%% Use the code for whatever you'd like                                %%%
%%%                                                                     %%%
%%% Demonstrates performing Total-Variation (TV) denoising and          %%%
%%% deblurring of a 1D input signal using the ADMM iterative scheme     %%%
%%% applied to a Generalized-Lasso problem.  More information on the    %%%
%%% approach can be found in [1].  The approach from [1] is modified    %%%
%%% slightly by using gradient-descent (Landweber iterations) to solve  %%%
%%% the first subproblem in place of a direct solver or iterative       %%%
%%% method such as conjugate gradient.  Matrix-free Landweber           %%%
%%% iterations are commonly used for large-scale linear inverse         %%%
%%% problems and this sample code allows the accuracy of the            %%%
%%% subproblem solves to be adjusted to see the effect on the final     %%%
%%% reconstructions.                                                    %%%
%%%                                                                     %%%
%%% [1] Boyd et al., Distributed Optimization and Statistical Learning  %%%
%%%     via the Alternating Direction Method of Multipliers, 2010       %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
close all;
clear all;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Problem Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
N           = 100;     % Signal sample count
lambda      = 0.2;     % Total-variation weight
sigma       = 9.0;     % PSF sigma for generating blurred input
noise       = 0.05;    % Gaussian-noise sigma
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADMM Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
rho         = 1.0;     % ADMM constraint weight 0 < rho <= 2
outer_iters = 100;     % Number of ADMM iterations
inner_iters = 2;       % Number of Landweber steps per ADMM iteration
relax       = 1.9;     % Under-relaxation factor for Landweber steps [0,2]
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct System PSF & Image Formation Model %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
off = -ceil(sigma*3):ceil(sigma*3);  % psf pixel offsets
psf = exp( -off.^2/sigma^2 );        % psf values for offsets
psf = psf/sum(psf);                  % normalize to 1.0
 
% generate psf matrix by setting diagonals based on 1D psf above
M = zeros( N, N );                   
for i=1:numel(psf),
    M = M + diag( psf(i)*ones(N-abs(off(i)),1), off(i) );
end
M = sparse( M );
 

% M=generate_PSF_matrix(sigma,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate noisy and blurred synthetic input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
input = zeros( N, 1 );                  % start with zeros
input(floor(N/4):ceil(3*N/4),:) = 1.0;  % make a center region at 1.0
% plot( input, 'b-+' );
% hold on;
% blur = M * input + noise*(randn(N,1));  % blur the input by the psf
% blur = M * input;  % blur the input by the psf
blur = input + noise*(randn(N,1));  % blur the input by the psf
% M=eye(size(M)); % TODO: 测试不加blur时，是否会出现smooth
% plot(blur,'.-g');
% hold off;
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct the difference matrix that computes image gradients %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
A = sparse( -diag( ones(N,1), 0 ) + diag( ones(N-1,1), 1 ) );
A(N,N-1)=1; % use a backward difference for the final point

% 测试自己写的带H的版本
lambda=0.3; % tvl2_h lambda=0.2, tvl1_h lambda=1
% xx=tvl1_total_variation_vec_H(blur,lambda,A,M);
xx=tvl2_total_variation_vec(blur,lambda,A);
figure;
hold on;
% plot( blur, 'r-o' );
plot( xx, 'b-o' );
axis tight;
axis([0 100 -0.2 1.2]);
set(gcf, 'Color', 'w');
hold off;
figure;
hold on;
plot( blur, 'r-o' );
% plot( xx, 'b-o' );
axis tight;
axis([0 100 -0.2 1.2]);
set(gcf, 'Color', 'w');
hold off;
return;

% 测试自己代码是否会有模糊的效果
% lambda=0.2;
% xx=tvl2_total_variation_vec_H(blur,lambda,A,M);
% x=tvl2_total_variation_vec(blur,lambda,A);


 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup ADMM and perform iterations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda      = 0.2;     % Total-variation weight
 
% compute the eigenvalues of the first subproblem
% system matrix and use the largest to define a
% Landweber step size, for large system a a power-
% iteration should be used instead.
tmp = eigs(M'*M + rho*A'*A);
step = relax/tmp(1);
 
% initialize the ADMM variables
% x = blur;             % intrinsic (sharp) signal
x = zeros( N, 1 );             % intrinsic (sharp) signal
z = zeros( N, 1 );    % splitting variable
u = zeros( N, 1 );    % scaled Lagrange multipliers
 
% define an anonymous shrinkage operator to implement
% the second sub-problem solve
shrink = @(kappa,x) max( abs(x)-kappa, 0 ).*sign(x);
 
% perform outer_iters outer iterations, plot the
% solver progress as the iterations proceed
for k=1:outer_iters,
    fprintf( 1, 'iteration %d of %d\n', k, outer_iters );
    % plot( x, 'b-+' );
    % drawnow;
  
    % define an anonymous function returning the gradient
    % of the first subproblem w.r.t. x, holding z and u
    % fixed, then perform gradient-descent (Landweber 
    % iterations).
    % gradF = @(x) M'*(M*x) - M'*blur + rho*A'*( A*x - z + u );
    % x = gradient_descent( gradF, x, step, inner_iters );
    xold=x;
    x=(M'*M+rho*A'*A)\(M'*blur+rho*A'*(z-u));
     
    % solve the second sub-problem using the anonymous
    % shrinkage operator
    zold=z;
    z = shrink( lambda/rho, A*x + u );
     
    % update the scaled Lagrange multipliers
    u = u + A*x - z;
    
    % convergence
    relres=norm(A*x-z);
    relchg=max(norm(z-zold),norm(x-xold));

    if relres<1e-4 && relchg<1e-4
        break;
    end
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup ADMM and perform iterations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
figure;
hold on;
% plot( blur, 'r-o' );
plot( x, 'b-o' );
plot( input, 'g-x' );
plot( xx, 'r-*' );
hold off;
