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
    
    if strcmp(options.nnm_mode, "real")
        manifold = rotationsfactory(n, 1);
        problem.M = manifold;
    else
        manifold = unitaryfactory(n, 1);
        problem.M = manifold;
    end
    
    if options.schur == true
        [Q_init, ~] = schur(A, options.nnm_mode);
    elseif strcmp(options.nnm_mode, "real")
        Q_init = randrot(n);
    else
        Q_init = randunitary(n);
    end
    

    problem.cost = @nnm_cost;
    problem.grad = @nnm_grad;
    problem.hess = @nnm_hess;  
    
    if options.testing_ == true
        checkgradient(problem);
        checkhessian(problem);
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
            if strcmp(options.nnm_mode, "real")
                store.S = PR(store.B, store.B);
            else
                store.S = P(store.B);
            end
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
        if strcmp(options.nnm_mode, "real")
            PCBU = PR(BU - UB, B);
        else
            PCBU = P(BU - UB);
        end
        H = 2 * skew(-UB' * S + BU * S' + comm(B', PCBU) - U * symm(store.almostgrad));
        % the formula from the paper is below 
        % H = 2 * skew(U * skew(comm(B', S)) + comm(comm(B, U)', S) + comm(B', P(comm(B, U))));
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

    function PW = PR(W, X)
        % projection "outside of quasidiagonal"
        % project W according to the branch of projection in X
        PW = zeros(size(W, 1));
        for i = 1:2:n-1
            X11 = X(i,i);
            X12 = X(i, i + 1);
            X21 = X(i + 1, i);
            X22 = X(i + 1,i + 1);
            if abs(X12 - X21) < abs(X11 - X22)
                PW(i,i) = W(i,i);
                PW(i + 1, i + 1) = W(i + 1, i + 1);
            else
                PW(i,i) = 0.5 * (W(i, i) + W(i + 1, i + 1));
                PW(i, i + 1) = 0.5 * (W(i, i + 1) - W(i + 1, i));
                PW(i + 1, i) = - PW(i, i + 1);
                PW(i + 1, i + 1) = PW(i, i);
            end
            if i == n - 2
                PW(n,n) = W(n,n);
            end
        end
        PW = W - PW;
    end
    
    % Below are functions that are currently not in use
    % function [cost, gradient] = nnm_costgrad(Q)
    %     B = Q' * A * Q;
    %     S = P(B);
    %     cost = norm(S, 'fro')^2;
    %     if nargout == 2
    %         G = 2 * (A * Q * S' + A' * Q * S);
    %         gradient = (Q' * G - G' * Q) / 2;
    %     end
    % end
    % 
    % function G = nnm_egrad(Q)
    %     B = Q' * A * Q;
    %     S = P(B);
    %     G = 2 * (A * Q * S' + A' * Q * S);
    % end
    % 
    % function eH = nnm_ehess(Q, U)
    %     Z = problem.M.tangent2ambient(Q, U);
    %     S = P(Q' * A * Q);
    %     dB = Z' * A * Q + Q' * A * Z;
    %     eH = A * Z * S' + A * Q * P(dB)' + A' * Z * S + A' * Q * P(dB);
    %     eH = 2 * eH;
    % end

end