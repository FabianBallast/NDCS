function [Aeq,beq, Af] = findCon(A,B, Tf, x0)
%FINDCON Summary of this function goes here
%   Detailed explanation goes here
n = length(A);
m = length(B(1,:));

Aeq = [A, B, -eye(n), zeros(n, (Tf-2)*(n+m) + m);
       zeros(n, n+m), A, B, -eye(n),  zeros(n, (Tf-3)*(n+m) + m);
       zeros(n, 2*(n+m)), A, B, -eye(n), zeros(n, (Tf-4)*(n+m) + m);
       zeros(n, 3*(n+m)), A, B, -eye(n), zeros(n, m);
       eye(n), zeros(n, (Tf-1)*(n+m) + m)];
Af = [zeros(n, 4*(n+m)), A, B];   
beq = [zeros((Tf-1)*n, 1); x0];
end

