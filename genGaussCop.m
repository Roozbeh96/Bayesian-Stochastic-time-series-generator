function U = genGaussCop(r, n)
rho = r;%2 * sin(r * pi/6); % Pearson correlation
P = toeplitz([1 rho]); % Correlation matrix
d = size(P, 1); % Dimension
% Generate sample
U = mvnrnd(zeros(d,1), P, n);
U = normcdf(U);
end