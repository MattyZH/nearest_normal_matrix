function U = random_unitary_haar(n)
    A = randn(n) + 1i*randn(n);
    [U, ~, V] = svd(A);       % U and V are unitary, singular values on diagonal of ~
    U = U * V';               % This combination is Haar-distributed unitary
end