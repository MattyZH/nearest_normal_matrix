function G = nnm_egrad(Q, A)
B = Q' * A * Q;
D = diag(diag(B));
S = B - D;
G = 2 * (A * Q * S' + A' * Q * S);
end
 