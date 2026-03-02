function  [A, dist_] = normal_matrix_noise(n, max_ratio)
% Creating a random normal matrix Atrue, then its noised version A,
% and returning A and norm(A - Atrue, 'fro')^2 as distance.
        Qt = randunitary(n);
        Dt = random_diagonal(n);
        Atrue = Qt * Dt * Qt';
        ratio = max_ratio * rand;
        A = add_noise(Atrue, ratio);
        dist_ = norm(A - Atrue, 'fro')^2;
end