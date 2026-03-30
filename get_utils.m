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
    utils.unitarybase2normalbase = @unitarybase2normalbase;
    utils.tangent_normal_projection = @tangent_normal_projection;
    utils.orthogonal_noise = @orthogonal_noise;
    utils.ortho_basis = @ortho_basis;
    utils.skew_hermitian_basis = @skew_hermitian_basis;


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

    function  [A, A_true, dist_] = normal_matrix_noise(n_, noise_ratio, noise_mode, is_noise_ortho)
        % Add noise to the diagonalization of A with norm 
        % noise_ratio * norm(A). 
        % noise_mode specifies if noise should be real or complex.
        % is_noise_ortho = true adds noise orthogonal to the manifold of
        % normal matrices in A.
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
        if is_noise_ortho
            Dt_noised = add_noise(Dt, noise_ratio, "ortho");
        else
            Dt_noised = add_noise(Dt, noise_ratio, noise_mode);
        end
        A = Qt * Dt_noised * Qt';
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

    function noise = orthogonal_noise(D)
        % For a diagonal matrix D, create a noise in a normal space of
        % the manifold of normal matrices in D.
        n = size(D, 1);
        basis = unitarybase2normalbase(D);
        basis = ortho_basis(basis); % checked: really orthonormal and tangent
        
        noise = randn(n) + 1i * randn(n);
        tangent_noise = tangent_normal_projection(noise, basis);
        noise = noise - tangent_noise;
    end

    function new_basis = ortho_basis(basis)
        % From a matrix-form basis to a matrix-form Frobenius-orthogonal basis.
        n = size(basis, 2);
        new_basis = zeros(size(basis)); 
        basis_vectorization = zeros(size(basis, 1), 2 * n * n);    
        for i = 1:size(basis, 1)
            X = squeeze(basis(i, :, :));
            rX = reshape(real(X), 1, n * n);
            iX = reshape(imag(X), 1, n * n);
            basis_vectorization(i, :) = [rX, iX];
        end
        ortho_basis_vectorization = transpose(orth(transpose(basis_vectorization)));
        for i = 1:size(ortho_basis_vectorization, 1)
            rY = reshape(ortho_basis_vectorization(i, 1:n*n), n, n);
            iY = reshape(ortho_basis_vectorization(i, n*n+1:end), n, n);
            new_basis(i, :, :) = rY + 1i * iY;
        end
    end

    function basis = skew_hermitian_basis(n)
        % Returns a basis of the space of n x n skew-hermitian matrices.
        % This is a basis of T_I U(n).
        basis = zeros(n * n - n, n, n);
    
        q = 0;
        sep = round(n * (n - 1) / 2);
    
        for i = 1:n-1
            for j = i+1:n
                q = q + 1;
                el1 = zeros(n,n);
                el1(i, j) = 1;
                el1(j, i) = -1;
                basis(q, :, :) = el1;
                
                el2 = zeros(n,n);
                el2(i, j) = 1i;
                el2(j, i) = 1i;
                basis(q + sep, :, :) = el2;
            end
        end
    end

    function projected_noise = tangent_normal_projection(noise, basis)
        % Projects the noise onto the subspace spanned by orthonormal
        % basis.
        rp_noise = zeros(size(noise));
        ip_noise = zeros(size(noise));
        for i = 1:size(basis, 1)
            X = squeeze(basis(i, :, :));
            rX = real(X);
            iX = imag(X);
            inner_product_ = real(trace(noise * X'));
            rp_noise = rp_noise + inner_product_ * rX;
            ip_noise = ip_noise + inner_product_ * iX;
        end
    
        projected_noise = rp_noise + 1i * ip_noise;
    end

    function normal_basis = unitarybase2normalbase(D)
        % Returns a basis of the tangent space of the normal matrix
        % manifold in D. D must be a diagonal matrix to work properly.
        % TODO add a warning if D is not diagonal.
        n = size(D, 1);
        u_basis = skew_hermitian_basis(n);
        normal_basis = zeros(n * (n + 1), n, n);
    
        for i = 1:n
            el = zeros(n);
            el(i, i) = D(i, i);
            normal_basis(i, :, :) = el;
    
            el = zeros(n);
            a = real(D(i,i));
            b = imag(D(i,i));
            el(i,i) = b - a * 1i;
            normal_basis(n + i, :, :) = el;
        end
    
        for j = 1:size(u_basis, 1)
            X = squeeze(u_basis(j, :, :));
            normal_basis(2 * n + j, :, :) = D * X + X' * D;
        end
    end
end