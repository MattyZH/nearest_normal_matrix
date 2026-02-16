function A_noised = add_noise(A, ratio)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
n = size(A);
noise = randn(n) + 1i * randn(n);
noise = noise / norm(noise, 'fro');
A_noised = A + ratio * noise * norm(A, 'fro');
end