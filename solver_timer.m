function [timing, timing_median, iterations] = solver_timer(solvers, options, test_data, distances, solver_up)
    
    verbosity = options.timer_verbosity;

    is_ = size(test_data, 1);
    batch_size = size(test_data{1}, 1);

    timing = zeros(size(solvers, 2), is_, batch_size);
    iterations = zeros(size(solvers, 2), is_, batch_size);
    for j = 1:numel(solvers)
        solver_ = solvers{j};
        if verbosity > 0
            fprintf('Started testing solver %s\n', func2str(solver_))
        end
        for i = 1:is_
            n = size(test_data{i}, 2);  
            if verbosity > 1
                fprintf('Started testing input size %d\n', n)
            end
            for k = 1:batch_size
                if verbosity > 2
                    fprintf('Started testing batch number %d out of %d. ', k, batch_size)
                end
                A = squeeze(test_data{i}(k, :, :));
                dist_ = distances{i}(k);

                % execution 
                [~, ~, iter_count, time_, cost_] = solver_up(A, solver_, options);
                
                if iter_count >= options.maxiter - 1
                    % too many iterations (9990+)
                    fprintf('Maximum iterations(%d) reached with solver %s (n=%d, batch=%d)\n', options.maxiter, func2str(solver_), n, k)
                end

                if cost_ > 1.1*dist_
                    % bad answer, found normal matrix is too far
                    fprintf('Cost exceeded threshold with solver %s (n=%d, batch=%d)\n', func2str(solver_), n, k);
                end

                iterations(j, i, k) = iter_count; 
                timing(j, i, k) = time_;
                if verbosity > 2
                    fprintf('Time = %f, %d iterations.\n', time_, iter_count)
                end
            end
        end
    end
    timing_median = median(timing, 3);

end