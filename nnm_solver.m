function [Qs, Qscost, iter_count, time_, cost_] = nnm_solver(A, solver, options)
%   
%   Finding the nearest normal matrix to A using the method 'solver'

n = size(A, 1);

manifold = unitaryfactory(n, 1);
problem.M = manifold;

Q_init = randunitary(n);
% Q_init = eye(n);

problem.cost = @(Q) nnm_cost(Q);
problem.costgrad = @(Q) nnm_costgrad(Q);
problem.egrad = @(Q) nnm_egrad(Q);
problem.ehess = @(Q, Z) nnm_ehess(Q, Z);
problem.hess = @(Q, Z) nnm_hess(Q, Z);

if options.testing_ == true
    checkhessian(problem);
end

[Qs, Qscost, info] = solver(problem, Q_init, options);
last_info = info(1, end);
iter_count = last_info.iter;
time_ = last_info.time;
cost_ = last_info.cost;

function cost = nnm_cost(Q)
    B = Q' * A * Q;
    S = P(B);
    cost = norm(S, 'fro')^2;
end


function [cost, gradient] = nnm_costgrad(Q)
    B = Q' * A * Q;
    S = P(B);
    cost = norm(S, 'fro')^2;
    if nargout == 2
        G = 2 * (A * Q * S' + A' * Q * S);
        gradient = (Q' * G - G' * Q) / 2;
    end
end

function G = nnm_egrad(Q)
    B = Q' * A * Q;
    S = P(B);
    G = 2 * (A * Q * S' + A' * Q * S);
end

function eH = nnm_ehess(Q, Z)
    Z = problem.M.tangent2ambient(Q, Z);
    S = P(Q' * A * Q);
    dB = Z' * A * Q + Q' * A * Z;
    eH = A * Z * S' + A * Q * P(dB)' + A' * Z * S + A' * Q * P(dB);
    eH = 2 * eH;
end

function H = nnm_hess(Q, Z)
    eH = nnm_ehess(Q, Z);
    G = nnm_egrad(Q);
    DG = eH - Q * Z * (Q' * G + G' * Q)/2;
    % % Z = Q * Z;
    % B = Q' * A * Q;
    % S = P(B);
    % % deltaB = comm(B, Z);
    % deltaB = Z' * Q' * A * Q + Q' * A * Q * Z;
    % deltaS = P(deltaB);
    % deltaQ = Q*Z;
    % DG = 2*(A * deltaQ * S' + A * Q * deltaS' + A' * deltaQ * S + A' * Q * deltaS);
    H = skew(Q' * DG);
end

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