function [Qs, Qscost, iter_count, time_, cost_] = nnm_solver(A, solver, options)
%NNM_SOLVER 
% Finding the 'normalization matrix' for A. If A1 is the nearest
% normal, then A1 = Qs' * D1 * Qs for some diagonal matrix D1.
%   
%
%   Inputs:
%       A        - complex-valued square matrix
%       SOLVER   - manopt solver for the problem (trustregions, arc, etc.)
%       OPTIONS  - struct of options for SOLVER
%
%   Outputs:
%       Qs        - 
%       Qscost    - norm(A - A1, 'fro')^2
%       iter_count- number of iterations executed by SOLVER 
%       time_     - time needed by SOLVER
%       cost_     - equals to Qscost
%
%   Example:
%       % 
%
%   See also
n = size(A, 1);

manifold = unitaryfactory(n, 1);
problem.M = manifold;

if options.schur == true
    [Q_init, ~] = schur(A);
else
    Q_init = randunitary(n);
end

problem.cost = @nnm_cost;
% problem.costgrad = @(Q) nnm_costgrad(Q);
problem.grad = @nnm_grad;
% problem.egrad = @(Q) nnm_egrad(Q);
% problem.ehess = @(Q, Z) nnm_ehess(Q, Z);
problem.hess = @nnm_hess;


if options.testing_ == true
    checkhessian(problem);
    checkgradient(problem)
end

[Qs, Qscost, info] = solver(problem, Q_init, options);
last_info = info(1, end);
iter_count = last_info.iter;
time_ = last_info.time;
cost_ = last_info.cost;

function store = prepare(Q, store)
    if ~isfield(store, 'B')
        store.B = Q' * A * Q;
    end
    if ~isfield(store, 'S')
        store.S = P(store.B);
    end
    if ~isfield(store, 'almostgrad')
        store.almostgrad = store.B * store.S' + store.B' * store.S;
    end

end

function [cost, store] = nnm_cost(Q, store)
    store = prepare(Q, store);
    cost = norm(store.S, 'fro')^2;
end

function [grad, store] = nnm_grad(Q, store)
    store = prepare(Q, store);
    grad = 2 * skew(store.almostgrad);
end

function [H, store] = nnm_hess(Q, U, store)
    % finally works
    store = prepare(Q, store);
    B = store.B;
    S = store.S;
    BU = B * U;
    UB = U * B;
    H = 2 * skew(-UB' * S + BU * S' + comm(B', P(BU - UB)) - U * symm(store.almostgrad));
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

function W = symm(W)
    W = 0.5 * (W + W');
end

% Below are functions that are currently not in use
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

function eH = nnm_ehess(Q, U)
    Z = problem.M.tangent2ambient(Q, U);
    S = P(Q' * A * Q);
    dB = Z' * A * Q + Q' * A * Z;
    eH = A * Z * S' + A * Q * P(dB)' + A' * Z * S + A' * Q * P(dB);
    eH = 2 * eH;
end

end