function [xopt] = findStates(A_large,B_large, A_T, B_T, x0, uopt)
%FINDSTATES Summary of this function goes here
%   Detailed explanation goes here
n = length(x0);
m = length(B_T(1, :));

Ae = [eye(n); A_large; A_T];
Be = [zeros(n, m); B_large; B_T];

xopt = Ae * x0 + Be * uopt;
end

