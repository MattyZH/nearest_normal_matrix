input_sizes = 16:16:128;% , 256, 512];
input_sizes = input_sizes(:);
batch_size = 20;
solvers = {@trustregions, @arc}; %, @conjugategradient};
max_ratio = 0.1;
[test_data, distances] = nnm_test_data_generator(input_sizes, ...
    batch_size, max_ratio);

options.verbosity = 0;
options.maxiter = 10000;
options.testing_ = false;
verbosity = 3;
[timing_rh, timing_average_rh, iterations_rh] = solver_timer(solvers, ...
    options, test_data, distances, verbosity);


plot(input_sizes, timing_average_rh, 'b-');

is_log = log(input_sizes(:));
tr_ta_log = log(tr_ta(:));
p = polyfit(is_log, tr_ta_log, 1);
alpha = p(1);
C = exp(p(2));

disp([alpha, C]);



