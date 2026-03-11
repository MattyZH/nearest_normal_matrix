[beta, gamma] = meshgrid(0:0.1:2 * pi, 0:0.05:pi);

A = randn(2) + 1i * randn(2);


Z = zeros(size(beta));
X = randn(1);
Y = randn(1);
for i=1:numel(Z)
    Z(i) = test_cost(1, gamma(i), A, beta(i) , Y(1));
end
% disp([beta, gamma])
surf(beta,gamma,Z)