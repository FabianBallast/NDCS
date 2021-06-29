function [A_large,B_large, A_T, B_T] = findLargeSystem(A,B,T)
%FINDLARGESYSTEM Return extended system matrices
%   Given a system x_k+1 = A * x_k + B * u_k, this returns:
%   [x_1;x_2;...;x_T-1] = A_large * x_0 + B_large * [u_0;u_1;...;u_T-2] 
% x_T = A_T * x_0 + B_T *  [u_0;u_1;...;u_T-1]
n = length(A);
m = length(B(1,:));
A_large = zeros(T*(n-1), n);
B_large = zeros(T*(n-1), T*m);

for i = 1:T-1
    A_large(1+n*(i-1):n*i, :) = A^i;
    for j = 1:i
        B_large(1+n*(i-1):n*i, 1+m*(j-1):m*j) = A^(i-j)*B;
    end
end

A_T = A^T;
B_T = zeros(n, T*m);
for j = 1:T
    B_T(:, 1+m*(j-1):m*j) = A^(T-j)*B;
end
end

