function [H, c] = findHessian(A_large,B_large, A_T, B_T, x0)
%FINDHESSIAN Summary of this function goes here
%   Detailed explanation goes here
n = length(x0);
m = length(B_T(1, :));

Ae = [eye(n); A_large]; %; A_T];
Be = [zeros(n, m); B_large]; %; B_T];

H = 2 * (Be.' * Be) + 2 * eye(length(Be(1, :)));
c = 2 * x0.' * Ae.' * Be;
end

