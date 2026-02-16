input_sizes = 128:16:256;% , 256, 512];
batch_size = 10;
solvers = {@trustregions}; %, @arc}; %, @conjugategradient};
options.verbosity = 0;
options.maxiter = 10000;
options.testing_ = false;
verbosity = 3;
[timing_new, timing_average_new, iterations_new] = solver_timer(solvers, options, input_sizes, batch_size, verbosity);

tr_ta = timing_average(:);
% tr_ta_new = timing_average_new(:);

plot(input_sizes, tr_ta, 'b-');

is_log = log(input_sizes(:));
tr_ta_log = log(tr_ta(:));
p = polyfit(is_log, tr_ta_log, 1);
alpha = p(1);
C = exp(p(2));

disp([alpha, C]);



