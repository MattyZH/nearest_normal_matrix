function [Q, cost_, info_] = nnm_solver(A, solver, options)
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
    
    utils = get_utils();

    comm = utils.commutator;
    skew = utils.skew;
    symm = utils.symm;
    P = utils.antidiag_projection;
    PR = utils.antiquasidiag_projection;

    n = size(A, 1);
    
    if strcmp(options.nnm_mode, "real")
        manifold = rotationsfactory(n, 1);
        problem.M = manifold;
    elseif strcmp(options.nnm_mode, "complex")
        manifold = unitaryfactory(n, 1);
        problem.M = manifold;
    else 
        error("Incorrect nnm_mode. Choose 'real' or 'complex'.")
    end
    
    if options.schur == true
        [Q_init, ~] = schur(A, options.nnm_mode);
    elseif strcmp(options.nnm_mode, "real")
        Q_init = randrot(n);
    elseif strcmp(options.nnm_mode, "complex")
        Q_init = randunitary(n);
    else 
    end
    

    problem.cost = @nnm_cost;
    problem.grad = @nnm_grad;
    problem.hess = @nnm_hess;  
    
    if options.testing_ == true
        checkgradient(problem);
        checkhessian(problem);
    end

    options.statsfun = @nnm_statsfun;
    
    [Q, cost_, info_] = solver(problem, Q_init, options);
    
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

    function stats = nnm_statsfun(~, Q, stats, store)
        if strcmp(options.nnm_mode, 'complex')
            current_approximation = Q * diag(diag(store.B)) * Q';
        else
            current_approximation = Q * (store.B - PR(store.B, store.B)) * Q';
        end
        stats.current_approx = current_approximation;
        stats.current_point = Q;
    end
end