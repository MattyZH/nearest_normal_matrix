input_sizes = 16:16:128;% , 256, 512];
input_sizes = input_sizes(:);
batch_size = 20;
solvers = {@arc}; %, @conjugategradient};
max_ratio = 0.01;
[test_data, distances] = nnm_test_data_generator(input_sizes, ...
    batch_size, max_ratio);

options.verbosity = 0;
options.maxiter = 10000;
options.testing_ = false;
options.schur = false;
verbosity = 3;
% warning('off', 'manopt:getHessian:approx');
[timing_random, timing_median_random, iterations_random] = solver_timer(solvers, ...
    options, test_data, distances, verbosity);
options.schur = true;
[timing_schur, timing_median_schur, iterations_schur] = solver_timer(solvers, ...
    options, test_data, distances, verbosity);


plot(input_sizes, timing_median_random(1, :), 'r-'); hold on;
plot(input_sizes, timing_median_schur(1, :), 'b-'); 
hold off;
legend('Random', 'Schur')


is_log = log(input_sizes(:));
tr_ta_log = log(timing_average_rh(1, :));
p = polyfit(is_log, tr_ta_log, 1);
alpha = p(1);
C = exp(p(2));

disp([alpha, C]);



