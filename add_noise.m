function A_noised = add_noise(A, ratio)
%   adding noise to A with Frobenius norm equal to 
%   RATIO * norm(A, 'fro'). The noise is uniformly distributed on a sphere.
%   Detailed explanation goes here
n = size(A);
noise = randn(n) + 1i * randn(n);
noise = noise / norm(noise, 'fro');
A_noised = A + ratio * noise * norm(A, 'fro');
end