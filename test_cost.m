function tc = test_cost(alpha, gamma, A, beta, delta)
Q = exp(-1i*alpha) * [exp(1i / 2  * (-delta - beta))*cos(gamma) -exp(1i / 2  * (delta - beta))*sin(gamma) ;
    exp(1i / 2  * (-delta + beta))*sin(gamma) exp(1i / 2  * (delta + beta))*cos(gamma)];
B = Q' * A * Q;
S = B - diag(diag(B));
tc = norm(B' * S - S * B', 'fro')^2;
end