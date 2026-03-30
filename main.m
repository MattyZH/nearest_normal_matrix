input_sizes = 128:1:128;% , 256, 512];
input_sizes = input_sizes(:);
batch_size = 10;
solver = @trustregions; %, @conjugategradient};
max_ratio = 0.03;

[test_data4, distances4] = nnm_test_data_generator(input_sizes, ...
    batch_size, max_ratio, "complex", false);
options.nnm_mode = "complex"; 
options.verbosity = 0;
options.maxiter = 10000;
options.testing_ = false;
options.schur = true;

options.timer_verbosity = 3;
test7_output = solver_timer(solver, options, test_data4, distances4, @nnm_solver);

options.nnm_mode = "real";
test4_output = solver_timer(solver, options, test_data, distances, @nnm_solver);

c_median = median(test_output.timing, 2);
r_median = median(test2_output.timing, 2);

plot(1:11, test_output.timing, 'r-'); hold on;
plot(1:11, test2_output.timing, 'b-'); 
hold off;
legend('Complex', 'Real')


is_log = log(input_sizes(1:end-1));
tr_ta_log = log(it_median(1:end-1));
p = polyfit(is_log, tr_ta_log, 1);
alpha = p(1);
C = exp(p(2));

disp([alpha, C]);

polyfit(input_sizes(1:end), it_mean(1:end), 1)



