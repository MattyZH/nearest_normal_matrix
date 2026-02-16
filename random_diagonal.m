function D = random_diagonal(n)
    d = randn(n,1) + 1i * randn(n,1);   
    D = diag(d);
end