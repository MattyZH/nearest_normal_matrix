function solver_output = solver_timer(solver, options, test_data, distances, solver_up)

    solver_output.noised_matrices = test_data{1, :};
    solver_output.true_matrices = test_data{2, :};
    solver_output.distances = distances;
    
    verbosity = options.timer_verbosity;

    n_input_sizes = size(test_data, 2);
    batch_size = size(test_data{1}, 1);
    
    solver_output.timing = zeros(n_input_sizes, batch_size);
    solver_output.iterations = zeros(n_input_sizes, batch_size);
    solver_output.info = cell(n_input_sizes, batch_size);
    solver_output.answer = cell(n_input_sizes, batch_size);
    solver_output.convergence = zeros(n_input_sizes, batch_size);
    solver_output.precise_convergence = zeros(n_input_sizes, batch_size);

    if verbosity > 0
        fprintf('Started testing solver %s\n', func2str(solver))
    end
    for i = 1:n_input_sizes
        n = size(test_data{1, i}, 2);  
        if verbosity > 1
            fprintf('Started testing input size %d\n', n)
        end
        for k = 1:batch_size
            if verbosity > 2
                fprintf('Started testing batch number %d out of %d. ', k, batch_size)
            end
            A = squeeze(test_data{1, i}(k, :, :));
            dist_ = distances{i}(k);

            [~, cost_, info_] = solver_up(A, solver, options);

            last_info = info_(end);
            n_iterations = last_info.iter;
            time_ = last_info.time;

            if verbosity > 2
                fprintf('Time = %f, %d iterations.\n', time_, n_iterations)
            end

            solver_output.timing(i, k) = time_;
            solver_output.iterations(i, k) = n_iterations;
            solver_output.info{i, k} = info_;
            solver_output.answer{i, k} = last_info.current_approx;

            if n_iterations >= options.maxiter - 1
                % algorithm did not converge
                solver_output.convergence(i, k) = 0;
                fprintf('Maximum iterations(%d) reached with solver %s (n=%d, batch=%d)\n', options.maxiter, func2str(solver), n, k)
            else
                solver_output.convergence(i, k) = 1;
            end

            if cost_ > 1.01*dist_
                % found normal matrix further than the true matrix.
                solver_output.precise_convergence(i, k) = 0;
                if solver_output.convergence(i, k)
                    fprintf('Cost exceeded threshold with solver %s (n=%d, batch=%d)\n', func2str(solver), n, k);
                end
            else
                solver_output.precise_convergence(i, k) = 1;
            end
        end
    end
end