function H = nnm_hess(Q, Sigma, A)
B = Q' * A * Q;
S = P(B);
H = Q * skew(comm(comm(B, Sigma), S) + comm(B, P(comm(B, Sigma))));
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