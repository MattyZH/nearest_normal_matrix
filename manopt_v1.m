% Generate random problem data.
n = 10;
A = random_unitary_haar(n);
% A = 11 * eye(n);

% Create the problem structure.
manifold = unitaryfactory(n, 1);
problem.M = manifold;

% Define the problem cost function and gradient
problem.costgrad = @(u) my_costgrad(A, u);

% Solve.
[x, xcost, info, options] = trustregions(problem, eye(n));

B = x' * A * x;
D = diag(diag(B));
disp(A)
A1 = x * D * x';
disp(A1)
disp(norm(A * A' - A' * A, 'fro')^2)
disp(norm(A1 * A1' - A1' * A1, 'fro')^2)
disp(norm(A - A1, 'fro')^2)


% disp(comm(x, x'))
