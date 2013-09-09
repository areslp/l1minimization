randn('seed', 0);
rand('seed',0);

m = 1500;       % amount of data
K = 5;        % number of blocks
partition = randi(50, [K 1]);

n = sum(partition); % number of features
p = 100/n;          % sparsity density

% generate block sparse solution vector
x = zeros(n,1);
start_ind = 1;
cum_part = cumsum(partition);
for i = 1:K,
    x(start_ind:cum_part(i)) = 0;
    if( rand() < p)
        % fill nonzeros
        x(start_ind:cum_part(i)) = randn(partition(i),1);
    end
    start_ind = cum_part(i)+1;
end