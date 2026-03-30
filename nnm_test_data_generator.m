function [test_data, distances] = nnm_test_data_generator(input_sizes, ...
    batch_sizes, noise_ratio, noise_mode, orthogonal_noise)

    utils = get_utils();

    normal_matrix_noise = utils.normal_matrix_noise; 

    is_ = size(input_sizes, 1);
    if size(batch_sizes, 1) == 1
        bs = batch_sizes;
        batch_sizes = bs * ones(is_, 1);
    end
    test_data = cell(2, is_);
    distances = cell(is_, 1);
    for i = 1:is_
        bs = batch_sizes(i);
        n = input_sizes(i);
        test_data{1, i} = zeros(bs, n, n);
        test_data{2, i} = zeros(bs, n, n);
        distances{i} = zeros(i, 1);
        for k = 1:bs
            [A, A_true, dist_] = normal_matrix_noise(n, noise_ratio, ...
                noise_mode, orthogonal_noise);
            distances{i}(k) = dist_;
            test_data{1, i}(k, :, :) = A;
            test_data{2, i}(k, :, :) = A_true;
        end
    end

end