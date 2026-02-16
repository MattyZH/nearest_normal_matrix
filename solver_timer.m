function [timing, timing_average, iterations] = solver_timer(solvers, options, input_sizes, batch_size, verbosity)
    timing = zeros(size(solvers, 2), size(input_sizes, 2), batch_size);
    iterations = zeros(size(solvers, 2), size(input_sizes, 2), batch_size);
    for j = 1:numel(solvers)
        solver_ = solvers{j};
        if verbosity > 0
            fprintf('Started testing solver %s\n', func2str(solver_))
        end
        for i = 1:numel(input_sizes)
            n = input_sizes(i);
            if verbosity > 1
                fprintf('Started testing input size %d\n', n)
            end
            for k = 1:batch_size
                if verbosity > 2
                    fprintf('Started testing batch number %d out of %d. ', k, batch_size)
                end
                % fprintf('Maximum iterations(%d) reached with solver %s (n=%d, batch=%d)\n', options.maxiter, func2str(solver), n, k)
                [A, dist_] = random_matrix(n);

                % execution 
                [Qs, Qscost, iter_count, time_, cost_] = nnm_solver(A, solver_, options);
                
                if iter_count >= options.maxiter
                    fprintf('Maximum iterations(%d) reached with solver %s (n=%d, batch=%d)\n', options.maxiter, func2str(solver_), n, k)
                end

                if cost_ > 2*dist_
                    fprintf('Cost exceeded threshold with solver %s (n=%d, batch=%d)\n', func2str(solver_), n, k);
                end

                iterations(j, i, k) = iter_count;
                timing(j, i, k) = time_; % Store the timing result
                if verbosity > 2
                    fprintf('Time = %f.\n', time_)
                end
            end
        end
    end
    timing_average = mean(timing, 3);
    
    function  [A, dist_] = random_matrix(n)
        Qt = randunitary(n);
        Dt = random_diagonal(n);
        Atrue = Qt * Dt * Qt';
        ratio = 0.05 * rand;
        A = add_noise(Atrue, ratio);
        dist_ = norm(A - Atrue, 'fro')^2;
    end
end