% Total variation denoising with random data

%% Generate problem data
rand('seed', 0);
randn('seed', 0);

n = 100;

% x0 = ones(n,1);
% for j = 1:3
    % idx = randsample(n,1);
    % k = randsample(1:10,1);
    % x0(ceil(idx/2):idx) = k*x0(ceil(idx/2):idx);
% end
% b = x0 + randn(n,1);

lambda = 1;

[y,t]=gen_test_signal;
b=y;
%% Solve problem
n=length(b);
e = ones(n,1);
D = spdiags([e -e], 0:1, n,n);
% x=tvl1_total_variation_vec(b',lambda,D); % lambda=2;
x=tvl2_total_variation_vec(b',lambda,D); % lambda=1;
% [x history] = total_variation_boyd(b', lambda, 1.0, 1.0);
hold on;
plot(t,x,'b-o');
hold off;

%% Reporting
% K = length(history.objval);                                                                                                        

% h = figure;
% plot(1:K, history.objval, 'k', 'MarkerSize', 10, 'LineWidth', 2); 
% ylabel('f(x^k) + g(z^k)'); xlabel('iter (k)');

% g = figure;
% subplot(2,1,1);                                                                                                                    
% semilogy(1:K, max(1e-8, history.r_norm), 'k', ...
    % 1:K, history.eps_pri, 'k--',  'LineWidth', 2); 
% ylabel('||r||_2'); 

% subplot(2,1,2);                                                                                                                    
% semilogy(1:K, max(1e-8, history.s_norm), 'k', ...
    % 1:K, history.eps_dual, 'k--', 'LineWidth', 2);   
% ylabel('||s||_2'); xlabel('iter (k)'); 
