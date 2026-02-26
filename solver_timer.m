function [timing, timing_average, iterations] = solver_timer(solvers, options, input_sizes, batch_size, verbosity)
%SOLVER_TIMER Times multiple solvers on varying input sizes.
%
%   [TIMING, TIMING_AVERAGE, ITERATIONS] = SOLVER_TIMER(SOLVERS, OPTIONS,
%       INPUT_SIZES, BATCH_SIZE, VERBOSITY) runs each solver in SOLVERS on
%       random problems of sizes specified in INPUT_SIZES, repeating each
%       size BATCH_SIZE times, and returns the timings, average timings, and
%       iteration counts.
%
%   Inputs:
%       SOLVERS      - Cell array of function handles to the solvers to be timed.
%       OPTIONS      - Struct containing optional parameters for the solvers.
%       INPUT_SIZES  - Array of input sizes to test the solvers on.
%       BATCH_SIZE   - Number of tests for each size to average over.
%       VERBOSITY    - Integer specifying the level of output verbosity (0 = silent, 1 = basic, etc.).
%
%   Outputs:
%       TIMING          - Array or reals containing raw timing data for each run.
%       TIMING_AVERAGE  - Array or reals containing average timings per solver and size.
%       ITERATIONS      - Array or integers containing iteration counts for each run.
%
%   Example:
%       solvers = {@solver1, @solver2};
%       options = struct('tol', 1e-6);
%       input_sizes = [10, 100, 1000];
%       batch_size = 5;
%       verbosity = 1;
%       [timing, avg, iters] = solver_timer(solvers, options, input_sizes, batch_size, verbosity);
%
%   See also TIC, TOC, TIMEIT.

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