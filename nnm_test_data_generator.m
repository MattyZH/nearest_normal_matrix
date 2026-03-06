function [test_data, distances] = nnm_test_data_generator(input_sizes, batch_sizes, max_ratio, mode)
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
            [A, dist_] = normal_matrix_noise(n, max_ratio, mode);
            distances{i}(k) = dist_;
            test_data{i}(k, :, :) = A;
        end
    end


    function  [A, dist_] = normal_matrix_noise(n_, max_ratio, mode)
        if strcmp(mode, "real")
            Qt = randrot(n_);
            Dt = random_quasidiagonal(n_);
        else
            disp(n_)
            Qt = randunitary(n_);
            Dt = random_diagonal(n_);
        end
        Atrue = Qt * Dt * Qt';
        ratio = max_ratio * rand;
        A = add_noise(Atrue, ratio, mode);
        dist_ = norm(A - Atrue, 'fro')^2;
        end


    function A_noised = add_noise(A, ratio, mode)
        sizeA = size(A);
        if strcmp(mode, "real")
            noise = randn(sizeA);
        else
            noise = randn(sizeA) + 1i * randn(sizeA);
        end
        noise = noise / norm(noise, 'fro');
        A_noised = A + ratio * noise * norm(A, 'fro');
    end

    function D = random_diagonal(n)
        d = randn(n, 1) + 1i * randn(n, 1);   
        D = diag(d);
    end
    function D = random_quasidiagonal(n)
        d = randn(n,1);  
        D = zeros(n);
        for j = 1:2:n-1
            D(j, j) = d(j); 
            D(j, j + 1) = d(j + 1);
            D(j + 1, j) = - d(j + 1);
            D(j + 1, j + 1) = d(j);
            if j == n - 2
                D(n, n) = d(n);
            end
        end
    end

end