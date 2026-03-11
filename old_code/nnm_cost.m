function cost = nnm_cost(Q, A)
B = Q' * A * Q;
D = diag(diag(B));
S = B - D;
cost = norm(S, 'fro')^2;
end