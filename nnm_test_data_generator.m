function [test_data, distances] = nnm_test_data_generator(input_sizes, batch_sizes, max_ratio)
    is_ = size(input_sizes, 1);
    if size(batch_sizes, 1) == 1
        bs = batch_sizes;
        batch_sizes = bs * ones(is_, 1);
    end
    test_data = cell(is_, 1);
    distances = cell(is_, 1);
    for i = 1:is_
        bs = batch_sizes(i);
        n = input_sizes(i);
        test_data{i} = zeros(i, n, n);
        distances{i} = zeros(i, 1);
        for k = 1:bs
            [A, dist_] = normal_matrix_noise(n, max_ratio);
            distances{i}(k) = dist_;
            test_data{i}(k, :, :) = A;
        end

    end

end