function [test_data, distances] = nnm_test_data_generator(input_sizes, ...
    batch_size, noise_ratio, noise_mode, is_ortho_noise)

    utils = get_utils();

    normal_matrix_noise = utils.normal_matrix_noise; 

    n_input_sizes = size(input_sizes, 1);
    test_data = cell(2, n_input_sizes);
    distances = cell(n_input_sizes, 1);
    for i = 1:n_input_sizes
        n = input_sizes(i);
        test_data{1, i} = zeros(batch_size, n, n);
        test_data{2, i} = zeros(batch_size, n, n);
        distances{i} = zeros(i, 1);
        for k = 1:batch_size
            [A, A_true, dist_] = normal_matrix_noise(n, noise_ratio, ...
                noise_mode, is_ortho_noise);
            distances{i}(k) = dist_;
            test_data{1, i}(k, :, :) = A;
            test_data{2, i}(k, :, :) = A_true;
        end
    end

end