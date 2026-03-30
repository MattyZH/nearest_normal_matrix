function utils = get_utils()
    % this file contains of utilitary functions that are used elsewhere.
    
    utils.commutator = @commutator;
    utils.skew = @skew;
    utils.antidiag_projection = @antidiag_projection;
    utils.antiquasidiag_projection = @antiquasidiag_projection;
    utils.symm = @symm;
    utils.random_diagonal = @random_diagonal;
    utils.random_quasidiagonal = @random_quasidiagonal;
    utils.add_noise = @add_noise;
    utils.normal_matrix_noise = @normal_matrix_noise;


    function W  = commutator(U, V)
        W = U * V - V * U;
    end
    
    function W = skew(W)
        W = 0.5 * (W - W');
    end
    
    function W = antidiag_projection(W)
        W = W - diag(diag(W));
    end
    
    function W = symm(W)
        W = 0.5 * (W + W');
    end

    function D = random_diagonal(n)
        % Create a random n x n complex-valued diagonal matrix.
        d = randn(n, 1) + 1i * randn(n, 1);   
        D = diag(d);
    end

    function D = random_quasidiagonal(n)
        % Create a random n x n NORMAL real-valued quasidiagonal matrix.
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

    function A_noised = add_noise(A, noise_ratio, noise_mode)
        % Create matrix A_noised with 
        % ||A_noised - A||_F = noise_ratio * ||A\||_F.
        % TODO noise_mode descriprion
        % TODO colliding eigvals in 'ortho' mode.
        sizeA = size(A);
        if strcmp(noise_mode, "ortho")
            if norm(A * A' - A' * A, 'fro')^2 > 1e-4
                warning("Matrix A should be normal for noise_mode = 'ortho'.")
            end
            noise = orthogonal_noise(A);
        elseif strcmp(noise_mode, "real")
            noise = randn(sizeA);
        elseif strcmp(noise_mode, "complex")
            noise = randn(sizeA) + 1i * randn(sizeA);
        else 
            error('Invalid noise mode specified. Choose "ortho", "real", or "complex".');
        end
        noise = noise / norm(noise, 'fro');
        A_noised = A + noise_ratio * noise * norm(A, 'fro');
    end

    function  [A, A_true, dist_] = normal_matrix_noise(n_, noise_ratio, noise_mode, orthogonal_noise)
        if strcmp(noise_mode, "real")
            Qt = randrot(n_);
            Dt = random_quasidiagonal(n_);
        elseif strcmp(noise_mode, 'complex')
            Qt = randunitary(n_);
            Dt = random_diagonal(n_);
        else
            error('Invalid noise mode specified. Choose "real" or "complex".')
        end
        A_true = Qt * Dt * Qt';
        % A_true = Dt;
        if orthogonal_noise
            Dt_noised = add_noise(Dt, noise_ratio, "ortho");
        else
            Dt_noised = add_noise(Dt, noise_ratio, noise_mode);
        end
        A = Qt * Dt_noised * Qt';
        % A = Dt_noised;
        % A = add_noise(A_true, noise_ration, mode_);
        dist_ = norm(A - A_true, 'fro')^2;
    end
    function PW = antiquasidiag_projection(W, X)
        % projection "outside of quasidiagonal"
        % project W according to the branch of projection in X
        % PW - W is a normal real-valued quasidiagonal matrix if X is
        % real-valued.
        PW = zeros(size(W, 1));
        for i = 1:2:n-1
            X11 = X(i,i);
            X12 = X(i, i + 1);
            X21 = X(i + 1, i);
            X22 = X(i + 1,i + 1);
            if abs(X12 - X21) < abs(X11 - X22)
                PW(i,i) = W(i,i);
                PW(i + 1, i + 1) = W(i + 1, i + 1);
            else
                PW(i,i) = 0.5 * (W(i, i) + W(i + 1, i + 1));
                PW(i, i + 1) = 0.5 * (W(i, i + 1) - W(i + 1, i));
                PW(i + 1, i) = - PW(i, i + 1);
                PW(i + 1, i + 1) = PW(i, i);
            end
            if i == n - 2
                PW(n,n) = W(n,n);
            end
        end
        PW = W - PW;
    end
end