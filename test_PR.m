function PW = test_PR(W, X)
        % projection "outside of quasidiagonal"
        % project W according to the branch of projection in X
        PW = zeros(size(W, 1));
        n = size(W, 1);
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
