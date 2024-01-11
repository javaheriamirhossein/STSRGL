function x = randl(m, n)
% Generate data from Laplace distribution

% m - number of rows
% n - number of columns
% x - data matrix with Laplacian distributed random samples 
%     with mean mu = 0 and std sigma = 1 (columnwise)


u1 = rand(m, n);
u2 = rand(m, n);

x = log(u1./u2);
x = bsxfun(@minus, x, mean(x));
x = bsxfun(@rdivide, x, std(x));

end