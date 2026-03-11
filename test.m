A = randn(6);% + 1i * randn(6);
disp(A)
options.testing_ = true;
options.tol = 1e-6;
options.verbosity = 0;
options.nnm_mode = 'real';
options.schur = true;
[Qs, Qscost, iter_count, time_, cost_] = nnm_solver(A, @arc, options);
B = Qs' * A * Qs;
D = B - test_PR(B, B);
A1 = Qs * D * Qs';
disp(A - A1)

