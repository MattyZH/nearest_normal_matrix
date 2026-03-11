function eH = nnm_ehess(Q, Z, A, problem)
Z = problem.M.tangent2ambient(Q, Z);
S = P(Q' * A * Q);
dB = Z' * A * Q + Q' * A * Z;
eH = A * Z * S' + A * Q * P(dB)' + A' * Z * S + A' * Q * P(dB);
eH = 2 * eH;
    function W  = comm(U, V)
        W = U * V - V * U;
    end
    function W = skew(W)
        W = 0.5 * (W - W');
    end
    function W = P(W)
        W = W - diag(diag(W));
    end
end