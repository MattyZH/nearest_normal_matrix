function [f, g] = my_costgrad(A,u)
    B = u' * A * u;
    D = diag(diag(B));
    f = norm(B, 'fro')^2 - norm(D, 'fro')^2;
    if nargout == 2
        g = +u * (comm(D', B) - ctranspose(comm(D', B)));
    end
end