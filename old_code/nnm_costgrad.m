function [cost, gradient] = nnm_costgrad(Q, A)
B = Q' * A * Q;
D = diag(diag(B));
S = B - D;
cost = norm(S, 'fro')^2;
if nargout == 2
    G = 2 * (A * Q * S' + A' * Q * S);
    gradient = (Q' * G - G' * Q) / 2;
end
 