function  convergence_graph(output_, data_)
% this function draw plots of convergence for the algorithm. 
% it expects output_ to be an output of solver_timer,
% and data_ to be an output of nnm_test_data_generator.
    distances_ = output_.distances;
    info_ = output_.info;
    true_matrices = data_(2, :);
    noised_matrices = data_(1, :);
    
    max_it = 0;
    
    for size_idx = 1:size(true_matrices, 2)
        for batch_n = 1:size(true_matrices{1, size_idx}, 1)
            n = size(true_matrices{1, size_idx}, 2);
    
            cur_info = info_{size_idx, batch_n};
            cur_it = size(cur_info, 2);
            cur_dist = distances_{size_idx}(batch_n);
            cur_normal = reshape(true_matrices{1, size_idx}(batch_n, :, :), n, n);
            cur_noised = reshape(noised_matrices{1, size_idx}(batch_n, :, :), n, n);
        
            scaled_distance_to_noised = zeros(cur_it, 1);
            scaled_distance_to_normal = zeros(cur_it, 1);
            scaled_distance_to_init = zeros(cur_it, 1);
        
            max_it = max(max_it, cur_it);
        
            for it = 1:cur_it
                cur_point = cur_info(it).current_approx;
                scaled_distance_to_noised(it) = norm(cur_point - cur_noised, "fro")^2 / cur_dist;
                scaled_distance_to_normal(it) = norm(cur_point - cur_normal, "fro")^2 / cur_dist;
                scaled_distance_to_init(it) = norm(cur_point - cur_info(1).current_approx, "fro")^2 / cur_dist;
            end
            plot(1:cur_it, scaled_distance_to_normal, 'b-'); hold on;
            plot(1:cur_it, scaled_distance_to_noised, 'g-'); hold on;
            % plot(1:cur_it, scaled_distance_to_init, 'y-'); hold on;
        end
    end
    plot(1:max_it, 0.5 * ones(1, max_it), 'r-');
    hold off;
    legend('scaled distance to normal', 'scaled distance to noised', 'scaled distance to init')
    xlabel('Iteration');
    ylabel('Scaled Distance');
    title('Comparison of Distances to Normal and Noised Matrices');
end