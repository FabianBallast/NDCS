%% clear old stuff, load data

clear; clc; close all;
aircraft

%% Find centralised solution
[AL1, BL1, AT1, BT1] = findLargeSystem(A1, B1, Tfinal);
[AL2, BL2, AT2, BT2] = findLargeSystem(A2, B2, Tfinal);
[AL3, BL3, AT3, BT3] = findLargeSystem(A3, B3, Tfinal);
[AL4, BL4, AT4, BT4] = findLargeSystem(A4, B4, Tfinal);

[H1, c1] = findHessian(AL1, BL1, AT1, BT1, x01);
[H2, c2] = findHessian(AL2, BL2, AT2, BT2, x02);
[H3, c3] = findHessian(AL3, BL3, AT3, BT3, x03);
[H4, c4] = findHessian(AL4, BL4, AT4, BT4, x04);

H = blkdiag(H1, H2, H3, H4);
c = [c1, c2, c3, c4];

ulim = umax / Tfinal * ones(length(H), 1);
zero_B = zeros(size(BT1));

Aeq = [BT1, -BT2, zero_B, zero_B;
       BT1, zero_B, -BT3, zero_B;
       BT1, zero_B, zero_B, -BT4];
   
beq = [AT2 * x02 - AT1 * x01;
       AT3 * x03 - AT1 * x01;
       AT4 * x04 - AT1 * x01];

options = optimoptions(@quadprog, 'Display', 'None');
[u_opt, fval,exitflag,output,lambda] = quadprog(H, c, [], [], Aeq, beq, -ulim, ulim, [], options);   
size_u = length(u_opt) / 4;

x_opt1 = findStates(AL1, BL1, AT1, BT1, x01, u_opt(1:size_u));
x_opt2 = findStates(AL2, BL2, AT2, BT2, x02, u_opt(size_u + 1:2*size_u));
x_opt3 = findStates(AL3, BL3, AT3, BT3, x03, u_opt(2*size_u +1:3*size_u));
x_opt4 = findStates(AL4, BL4, AT4, BT4, x04, u_opt(3*size_u +1:4*size_u));

n = length(A1);

figure(100);
plot(0:4, u_opt(1:2:size_u), '-o'); hold on;
plot(0:4, u_opt(2:2:size_u), '-x');  hold off;
legend('u_1', 'u_2');
grid on;
title('Optimal input airplane 1.');
xlabel('Timestep t');
ylabel('Input');

figure(101);
plot(0:4, u_opt(size_u + 1:2:2*size_u), '-o'); hold on;
plot(0:4, u_opt(size_u + 2:2:2*size_u), '-x');  hold off;
legend('u_1', 'u_2');
grid on;
title('Optimal input airplane 2.');
xlabel('Timestep t');
ylabel('Input');


figure(102);
plot(0:4, u_opt(2*size_u +1:2:3*size_u), '-o'); hold on;
plot(0:4, u_opt(2*size_u +2:2:3*size_u), '-x');  hold off;
legend('u_1', 'u_2');
grid on;
title('Optimal input airplane 3.');
xlabel('Timestep t');
ylabel('Input');


figure(103);
plot(0:4, u_opt(3*size_u +1:2:4*size_u), '-o'); hold on;
plot(0:4, u_opt(3*size_u +2:2:4*size_u), '-x');  hold off;
legend('u_1', 'u_2');
grid on;
title('Optimal input airplane 4.');
xlabel('Timestep t');
ylabel('Input');


for i = 1:n
    figure(i);
    plot(0:5, x_opt1(i:n:end), '-o'); hold on;
    plot(0:5, x_opt2(i:n:end), '-x'); 
    plot(0:5, x_opt3(i:n:end), '-d'); 
    plot(0:5, x_opt4(i:n:end), '-p'); hold off;
    legend('x_1', 'x_2', 'x_3', 'x_4');
    grid on;
    title(['State ', num2str(i)]);
    xlabel('Timestep t');
    ylabel('Position');
end
    
%% Question a dual
theta0 = [0;0;0;0];
mu_0 = [0;0;0;0];
alpha_0 = 0.005;
iterations = 15000;

theta_k = zeros(4, iterations+1, 4);
theta_k(:, 1, 1) = theta0;
theta_k(:, 1, 2) = theta0;
theta_k(:, 1, 3) = theta0;
theta_k(:, 1, 4) = theta0;
ulim = umax / Tfinal * ones(size_u, 1);

mu_k = zeros(4, iterations+1, 3);
mu_k(:, 1, 1) = mu_0;
mu_k(:, 1, 2) = mu_0;
mu_k(:, 1, 3) = mu_0;

c_ex1a = c1 + (mu_k(:, 1, 1)+mu_k(:, 1, 2)+mu_k(:, 1, 3)).' * BT1;
c_ex2 = c2 - mu_k(:, 1, 1).' * BT2;
c_ex3 = c3 - mu_k(:, 1, 2).' * BT3;
c_ex4 = c4 - mu_k(:, 1, 3).' * BT4;

for k = 1:iterations
    if (k > 8000)
        alpha_k=alpha_0/4;
    else
        alpha_k=alpha_0/1;
    end
    u_opt1 = quadprog(H1, c_ex1a, [], [], [], [], -ulim, ulim, [], options); 
    u_opt2 = quadprog(H2, c_ex2, [], [], [], [], -ulim, ulim, [], options); 
    u_opt3 = quadprog(H3, c_ex3, [], [], [], [], -ulim, ulim, [], options); 
    u_opt4 = quadprog(H4, c_ex4, [], [], [], [], -ulim, ulim, [], options); 
    
    theta_k(:, k+1, 1) = AT1 * x01 + BT1 * u_opt1;
    theta_k(:, k+1, 2) = AT2 * x02 + BT2 * u_opt2;
    theta_k(:, k+1, 3) = AT3 * x03 + BT3 * u_opt3;
    theta_k(:, k+1, 4) = AT4 * x04 + BT4 * u_opt4;

    mu_k(:, k+1, 1) = mu_k(:, k, 1) + alpha_k * (theta_k(:, k+1, 1)-theta_k(:, k+1, 2));
    mu_k(:, k+1, 2) = mu_k(:, k, 2) + alpha_k * (theta_k(:, k+1, 1)-theta_k(:, k+1, 3));
    mu_k(:, k+1, 3) = mu_k(:, k, 3) + alpha_k * (theta_k(:, k+1, 1)-theta_k(:, k+1, 4));
    
    c_ex1a = c1 + (mu_k(:, k+1, 1)+mu_k(:, k+1, 2)+mu_k(:, k+1, 3)).' * BT1;
    c_ex2 = c2 - mu_k(:, k+1, 1).' * BT2;
    c_ex3 = c3 - mu_k(:, k+1, 2).' * BT3;
    c_ex4 = c4 - mu_k(:, k+1, 3).' * BT4;
end

theta_opt = x_opt1(end-3:end);
error1 = theta_k(:, :, 1) - theta_opt;
rel_err1 = vecnorm(error1) / norm(theta_opt);
error2 = theta_k(:, :, 2) - theta_opt;
rel_err2 = vecnorm(error2) / norm(theta_opt);
error3 = theta_k(:, :, 3) - theta_opt;
rel_err3 = vecnorm(error3) / norm(theta_opt);
error4 = theta_k(:, :, 4) - theta_opt;
rel_err4 = vecnorm(error4) / norm(theta_opt);

rel_err = [rel_err1;
           rel_err2;
           rel_err3;
           rel_err4];

figure(5);
semilogy(0:iterations, rel_err(4, :));
title("Convergence of local estimate airplane 4.");
grid on;
xlabel('Iteration k');
ylabel('$$\frac{||\theta_4(k)-\theta_4^*||}{||\theta_4^*||}$$', 'Interpreter', 'Latex');

mus_opt = lambda.eqlin;
mu1_opt = mus_opt(1:4);
mu2_opt = mus_opt(5:8);
mu3_opt = mus_opt(9:12);

error1 = mu_k(:, :, 1) - mu1_opt;
rel_err_mu1 = vecnorm(error1) / norm(mu1_opt);
error2 = mu_k(:, :, 2) - mu2_opt;
rel_err_mu2 = vecnorm(error2) / norm(mu2_opt);
error3 = mu_k(:, :, 3) - mu3_opt;
rel_err_mu3 = vecnorm(error3) / norm(mu3_opt);

rel_err_mu = [rel_err_mu1;
           rel_err_mu2;
           rel_err_mu3];

figure(6);
semilogy(0:iterations, rel_err_mu(1, :)); hold on;
semilogy(0:iterations, rel_err_mu(2, :)); 
semilogy(0:iterations, rel_err_mu(3, :)); hold off;
title("Convergence of dual variables.");
legend('\mu_1', '\mu_2', '\mu_3');
grid on;
xlabel('Iteration k');
ylabel('$$\frac{||\mu_i(k)-\mu_i^*||}{||\mu_i^*||}$$', 'Interpreter', 'Latex');

xopt1 = findStates(AL1, BL1, AT1, BT1, x01, u_opt1);
xopt2 = findStates(AL2, BL2, AT2, BT2, x02, u_opt2);
xopt3 = findStates(AL3, BL3, AT3, BT3, x03, u_opt3);
xopt4 = findStates(AL4, BL4, AT4, BT4, x04, u_opt4);

n = length(A1);

for i = 1:n
    figure(6+i);
    plot(0:5, xopt1(i:n:end), '-o'); hold on;
    plot(0:5, xopt2(i:n:end), '-x'); 
    plot(0:5, xopt3(i:n:end), '-d'); 
    plot(0:5, xopt4(i:n:end), '-p'); hold off;
    legend('x_1', 'x_2', 'x_3', 'x_4');
    grid on;
    title(['State ', num2str(i)]);
    xlabel('Timestep t');
    ylabel('Position');
end

%% Question b primal
% Constants
alpha = [0.048, 0.044, 0.04, 0.03, 0.02, 0.01];
iterations = 100;
theta_k = zeros(4, iterations+1, length(alpha));
ulim = umax / Tfinal * ones(size_u, 1);
for k = 1:iterations
    for h = 1:length(alpha)
        [u_opt1, fval,exitflag,output,lambda1] = quadprog(H1, c1, [], [], BT1, theta_k(:, k, h) - AT1*x01, -ulim, ulim, [], options); 
        [u_opt2, fval,exitflag,output,lambda2] = quadprog(H2, c2, [], [], BT2, theta_k(:, k, h) - AT2*x02, -ulim, ulim, [], options); 
        [u_opt3, fval,exitflag,output,lambda3] = quadprog(H3, c3, [], [], BT3, theta_k(:, k, h) - AT3*x03, -ulim, ulim, [], options); 
        [u_opt4, fval,exitflag,output,lambda4] = quadprog(H4, c4, [], [], BT4, theta_k(:, k, h) - AT4*x04, -ulim, ulim, [], options); 

        theta_k(:, k+1, h) = theta_k(:, k, h) + alpha(h) * (lambda1.eqlin + lambda2.eqlin + lambda3.eqlin + lambda4.eqlin);
    end
end

theta_opt = x_opt1(end-3:end);
rel_err_con = zeros(length(alpha), iterations+1);
for h = 1:length(alpha)
    error = theta_k(:, :, h) - theta_opt;
    rel_err_con(h, :) = vecnorm(error) / norm(theta_opt);
end

figure(11);
semilogy(0:iterations, rel_err_con(1, :)); hold on;
semilogy(0:iterations, rel_err_con(2, :)); 
semilogy(0:iterations, rel_err_con(3, :)); 
semilogy(0:iterations, rel_err_con(4, :)); 
semilogy(0:iterations, rel_err_con(5, :));
semilogy(0:iterations, rel_err_con(6, :));hold off;
grid on;
legend('\alpha=0.048', '\alpha=0.044', '\alpha=0.04', '\alpha=0.03', '\alpha=0.02', '\alpha=0.01');
title("Convergence of target position for different constant step sizes.");
xlabel('Iteration k');
ylabel('$$\frac{||\theta(k)-\theta^*||}{||\theta^*||}$$', 'Interpreter', 'Latex');

% Variable, /k
alpha = [0.13, 0.1, 0.05, 0.01];
iterations = 100;
theta_k = zeros(4, iterations+1, length(alpha));
ulim = umax / Tfinal * ones(size_u, 1);
for k = 1:iterations
    for h = 1:length(alpha)
        alpha_k = alpha(h) / k;
        [u_opt1, fval,exitflag,output,lambda1] = quadprog(H1, c1, [], [], BT1, theta_k(:, k, h) - AT1*x01, -ulim, ulim, [], options); 
        [u_opt2, fval,exitflag,output,lambda2] = quadprog(H2, c2, [], [], BT2, theta_k(:, k, h) - AT2*x02, -ulim, ulim, [], options); 
        [u_opt3, fval,exitflag,output,lambda3] = quadprog(H3, c3, [], [], BT3, theta_k(:, k, h) - AT3*x03, -ulim, ulim, [], options); 
        [u_opt4, fval,exitflag,output,lambda4] = quadprog(H4, c4, [], [], BT4, theta_k(:, k, h) - AT4*x04, -ulim, ulim, [], options); 

        theta_k(:, k+1, h) = theta_k(:, k, h) + alpha_k * (lambda1.eqlin + lambda2.eqlin + lambda3.eqlin + lambda4.eqlin);
    end
end

theta_opt = x_opt1(end-3:end);
rel_err = zeros(length(alpha), iterations+1);
for h = 1:length(alpha)
    error = theta_k(:, :, h) - theta_opt;
    rel_err(h, :) = vecnorm(error) / norm(theta_opt);
end

figure(12);
semilogy(0:iterations, rel_err(1, :)); hold on;
semilogy(0:iterations, rel_err(2, :)); 
semilogy(0:iterations, rel_err(3, :)); 
semilogy(0:iterations, rel_err(4, :)); hold off;
grid on;
legend('\alpha_0=0.13', '\alpha_0=0.1', '\alpha_0=0.05', '\alpha_0=0.01');
title("Convergence of target position for different variable step sizes.");
xlabel('Iteration k');
ylabel('$$\frac{||\theta(k)-\theta^*||}{||\theta^*||}$$', 'Interpreter', 'Latex');


% Variable, /root
alpha = [0.13, 0.10, 0.085, 0.075];
iterations = 100;
theta_k = zeros(4, iterations+1, length(alpha));
ulim = umax / Tfinal * ones(size_u, 1);
for k = 1:iterations
    for h = 1:length(alpha)
        alpha_k = alpha(h) / nthroot(k, h);
        [u_opt1, fval,exitflag,output,lambda1] = quadprog(H1, c1, [], [], BT1, theta_k(:, k, h) - AT1*x01, -ulim, ulim, [], options); 
        [u_opt2, fval,exitflag,output,lambda2] = quadprog(H2, c2, [], [], BT2, theta_k(:, k, h) - AT2*x02, -ulim, ulim, [], options); 
        [u_opt3, fval,exitflag,output,lambda3] = quadprog(H3, c3, [], [], BT3, theta_k(:, k, h) - AT3*x03, -ulim, ulim, [], options); 
        [u_opt4, fval,exitflag,output,lambda4] = quadprog(H4, c4, [], [], BT4, theta_k(:, k, h) - AT4*x04, -ulim, ulim, [], options); 

        theta_k(:, k+1, h) = theta_k(:, k, h) + alpha_k * (lambda1.eqlin + lambda2.eqlin + lambda3.eqlin + lambda4.eqlin);
    end
end

theta_opt = x_opt1(end-3:end);
rel_err = zeros(length(alpha), iterations+1);
for h = 1:length(alpha)
    error = theta_k(:, :, h) - theta_opt;
    rel_err(h, :) = vecnorm(error) / norm(theta_opt);
end

figure(13);
semilogy(0:iterations, rel_err(1, :)); hold on;
semilogy(0:iterations, rel_err(2, :)); 
semilogy(0:iterations, rel_err(3, :)); 
semilogy(0:iterations, rel_err(4, :)); hold off;
grid on;
legend('$$\alpha(k)=\frac{0.13}{k}$$', '$$\alpha(k)=\frac{0.1}{k^{1/2}}$$', '$$\alpha(k)=\frac{0.085}{k^{1/3}}$$', '$$\alpha(k)=\frac{0.075}{k^{1/4}}$$', 'Interpreter', 'Latex');
title("Convergence of target position for different variable step sizes.");
xlabel('Iteration k');
ylabel('$$\frac{||\theta(k)-\theta^*||}{||\theta^*||}$$', 'Interpreter', 'Latex');

%% Question c Nesterov
theta0 = [0;0;0;0];
alpha = [  0.025 0.0275 0.03 0.0325  ];
iterations = 100;
options = optimoptions(@quadprog, 'Display', 'None');
theta_k = zeros(4, iterations+2, length(alpha));
lambdak = zeros(length(alpha), 1);
lambdap = zeros(length(alpha), 1);
for h = 1:length(alpha)
    theta_k(:, 1, h) = theta0;
end

ulim = umax / Tfinal * ones(size_u, 1);

for k = 1:iterations
    for h = 1:length(alpha) 
        lambdak(h) = (1+sqrt(1+4*lambdap(h)^2))/2;
        betak = (lambdap(h)-1)/lambdak(h);
        yk = theta_k(:, k+1, h) + betak*(theta_k(:, k+1, h)-theta_k(:, k, h));
        [u_opt1, fval,exitflag,output,lambda1] = quadprog(H1, c1, [], [], BT1, yk - AT1*x01, -ulim, ulim, [], options); 
        [u_opt2, fval,exitflag,output,lambda2] = quadprog(H2, c2, [], [], BT2, yk  - AT2*x02, -ulim, ulim, [], options); 
        [u_opt3, fval,exitflag,output,lambda3] = quadprog(H3, c3, [], [], BT3, yk  - AT3*x03, -ulim, ulim, [], options); 
        [u_opt4, fval,exitflag,output,lambda4] = quadprog(H4, c4, [], [], BT4, yk  - AT4*x04, -ulim, ulim, [], options); 
        theta_k(:, k+2, h) = yk + alpha(h) * (lambda1.eqlin + lambda2.eqlin + lambda3.eqlin + lambda4.eqlin);
        
        if(norm(theta_k(:, k+2, h)-theta_opt) > norm(theta_k(:, k+1, h)-theta_opt))
            theta_k(:, k+2, h) = theta_k(:, k+1, h);
            lambdak(h) = 0;
        end
        lambdap(h) = lambdak(h);
    end
    
end

theta_opt = x_opt1(end-3:end);
rel_err = zeros(length(alpha), iterations+1);
for h = 1:length(alpha)
    error = theta_k(:, 2:end, h) - theta_opt;
    rel_err(h, :) = vecnorm(error) / norm(theta_opt);
end

figure(14); 
for h = 1:length(alpha) 
    semilogy(0:iterations, rel_err(h, :)); hold on;
end
hold off;
grid on;
legend('\alpha=0.025', '\alpha=0.0.275', '\alpha=0.03', '\alpha=0.0325', '\alpha=0.035');
title("Convergence of target position with Nesterov method.");
xlabel('Iteration k');
ylabel('$$\frac{||\theta(k)-\theta^*||}{||\theta^*||}$$', 'Interpreter', 'Latex');

%% Question d
W = [0.75 0.25 0 0;
     0.25 0.50 0.25 0;
     0 0.25 0.5 0.25;
     0 0 0.25 0.75];

theta0 = [0;0;0;0];
alpha_0 = 0.03*4;
iterations = 100;
phi = [0 1 10 30 100];

W_phi = zeros(4, 4, length(phi));
for i = 1:length(phi)
    W_phi(:, :, i) = W^phi(i);
end

theta_k = zeros(4, iterations+1, 4, length(phi));
for i = 1:4
    for j = 1:length(phi)
        theta_k(:, 1, i, j) = theta0;
    end
end

ulim = umax / Tfinal * ones(size_u, 1);
for k = 1:iterations
    for h = 1:length(phi)
        [u_opt1, fval,exitflag,output,lambda1] = quadprog(H1, c1, [], [], BT1, theta_k(:, k, 1, h) - AT1*x01, -ulim, ulim, [], options); 
        [u_opt2, fval,exitflag,output,lambda2] = quadprog(H2, c2, [], [], BT2, theta_k(:, k, 2, h) - AT2*x02, -ulim, ulim, [], options); 
        [u_opt3, fval,exitflag,output,lambda3] = quadprog(H3, c3, [], [], BT3, theta_k(:, k, 3, h) - AT3*x03, -ulim, ulim, [], options); 
        [u_opt4, fval,exitflag,output,lambda4] = quadprog(H4, c4, [], [], BT4, theta_k(:, k, 4, h) - AT4*x04, -ulim, ulim, [], options); 

        lambdas = [lambda1.eqlin, lambda2.eqlin, lambda3.eqlin, lambda4.eqlin];

        for i = 1:4
            sum = [0;0;0;0];
            for j = 1:4
                sum = sum + W_phi(i, j, h) * (theta_k(:, k, j, h) + alpha_0 * lambdas(:, j));
            end
            theta_k(:, k+1, i, h) = sum;
        end
    end
end

theta_opt = x_opt1(end-3:end);
rel_err = zeros(length(phi), iterations+1, 4);
for j = 1:4
    for h = 1:length(phi)
        error = theta_k(:, :, j, h) - theta_opt;
        rel_err(h, :, j) = vecnorm(error) / norm(theta_opt);
    end


    figure(16+j); 
    for h = 1:length(phi) 
        semilogy(0:iterations, rel_err(h, :, j)); hold on;
    end
    semilogy(0:iterations, rel_err_con(4, :), '--');
    hold off;
    grid on;
    legend('\phi=0', '\phi=1', '\phi=10', '\phi=30', '\phi=100', 'Standard subgradient');
    title(['Convergence of target position of airplane ', num2str(j), ' with Consensus approach.']);
    xlabel('Iteration k'); 
    ylabel('$$\frac{||\theta(k)-\theta^*||}{||\theta^*||}$$', 'Interpreter', 'Latex');
end
%% Question e with fmincon

fval= @(x) 0.5*x.'*H*x + c*x;
lb = -(umax) * ones(length(H), 1);
ub =  (umax) * ones(length(H), 1);
options = optimoptions('fmincon', 'MaxFunctionEvaluations', 30000, 'StepTolerance', 1e-14);
x = fmincon(fval, u_opt, [], [], Aeq, beq, lb, ub, @nonlcon, options);

xopt1 = findStates(AL1, BL1, AT1, BT1, x01, x(1:10));
xopt2 = findStates(AL2, BL2, AT2, BT2, x02, x(11:20));
xopt3 = findStates(AL3, BL3, AT3, BT3, x03, x(21:30));
xopt4 = findStates(AL4, BL4, AT4, BT4, x04, x(31:40));

n = length(A1);

for i = 1:n
    figure(20+i);
    plot(0:5, xopt1(i:n:end), '-o'); hold on;
    plot(0:5, xopt2(i:n:end), '-x'); 
    plot(0:5, xopt3(i:n:end), '-d'); 
    plot(0:5, xopt4(i:n:end), '-p'); hold off;
    legend('x_1', 'x_2', 'x_3', 'x_4');
    grid on;
    title(['State ', num2str(i)]);
    xlabel('Timestep t');
    ylabel('Position');
end
%% Question e with Gurobi
% Centralised solution 

[Aeq1, Beq1, Af1] = findCon(A1, B1, Tfinal, x01);
[Aeq2, Beq2, Af2] = findCon(A2, B2, Tfinal, x02);
[Aeq3, Beq3, Af3] = findCon(A3, B3, Tfinal, x03);
[Aeq4, Beq4, Af4] = findCon(A4, B4, Tfinal, x04);

size_H = length(Aeq1(1, :));
H_t= eye(4*size_H);

zero_A = zeros(size(Aeq1));
zero_Af=  zeros(size(Af1));
Aeq_t = [Aeq1, zero_A, zero_A, zero_A;
        zero_A, Aeq2, zero_A, zero_A;
        zero_A, zero_A, Aeq3, zero_A;
        zero_A, zero_A, zero_A, Aeq4;
        Af1, -Af2, zero_Af, zero_Af;
        Af1, zero_Af, -Af3, zero_Af;
        Af1, zero_Af, zero_Af, -Af4];

beq_t = [Beq1; Beq2; Beq3; Beq4; zeros(12, 1)];
model.Q = sparse(H_t);
model.A = sparse(Aeq_t);
model.rhs = beq_t.';
model.sense = '=';

params.NumericFocus = 3;
params.outputflag = 0;
params.QCPdual = 1;

zero_H = zeros(size(H_t)/4);
ones_H = blkdiag(zeros(4), eye(2));
ones_H = blkdiag(ones_H, ones_H, ones_H, ones_H, ones_H);
model.quadcon(1).Qc = sparse(blkdiag(ones_H, zero_H, zero_H, zero_H));
model.quadcon(2).Qc = sparse(blkdiag(zero_H, ones_H, zero_H, zero_H));
model.quadcon(3).Qc = sparse(blkdiag(zero_H, zero_H, ones_H, zero_H));
model.quadcon(4).Qc = sparse(blkdiag(zero_H, zero_H, zero_H, ones_H));

model.quadcon(1).q = sparse(zeros(length(H_t), 1));
model.quadcon(2).q = sparse(zeros(length(H_t), 1));
model.quadcon(3).q = sparse(zeros(length(H_t), 1));
model.quadcon(4).q = sparse(zeros(length(H_t), 1));

model.quadcon(1).rhs = umax^2;
model.quadcon(2).rhs = umax^2;
model.quadcon(3).rhs = umax^2;
model.quadcon(4).rhs = umax^2;


ub_s = [400; 400; 400; 400; umax; (umax)];
ub_L = [ub_s; ub_s; ub_s; ub_s; ub_s];
ub_L = [ub_L; ub_L; ub_L; ub_L];

model.lb = -ub_L;
model.ub = ub_L;

result = gurobi(model, params);

x1 = result.x(1:30);
x2 = result.x(31:60);
x3 = result.x(61:90);
x4 = result.x(91:120);
u_pos = [5, 6, 11, 12, 17, 18, 23, 24, 29, 30];
xopt1 = findStates(AL1, BL1, AT1, BT1, x01, x1(u_pos));
xopt2 = findStates(AL2, BL2, AT2, BT2, x02, x2(u_pos));
xopt3 = findStates(AL3, BL3, AT3, BT3, x03, x3(u_pos));
xopt4 = findStates(AL4, BL4, AT4, BT4, x04, x4(u_pos));

n = length(A1);

for i = 1:n
    figure(24+i);
    plot(0:5, xopt1(i:n:end), '-o'); hold on;
    plot(0:5, xopt2(i:n:end), '-x'); 
    plot(0:5, xopt3(i:n:end), '-d'); 
    plot(0:5, xopt4(i:n:end), '-p'); hold off;
    legend('x_1', 'x_2', 'x_3', 'x_4');
    grid on;
    title(['State ', num2str(i)]);
    xlabel('Timestep t');
    ylabel('Position');
end

%% Question E dual 
iterations = 200;
theta0 = [0;0;0;0];
theta_k = zeros(4, iterations+1);
theta_k(:, 1) = theta0;
alpha_k = 0.044;
ub_s = [400; 400; 400; 400; umax; umax];
ub_L = [ub_s; ub_s; ub_s; ub_s; ub_s];

model1.Q = sparse(eye(size_H));
model1.A = sparse([Aeq1; Af1]);
model1.rhs = [Beq1; theta0];
model1.sense = '=';
model1.quadcon(1).Qc = sparse(ones_H);
model1.quadcon(1).q = sparse(zeros(length(H_t)/4, 1));
model1.quadcon(1).rhs = umax^2;
model1.lb = -ub_L;
model1.ub = ub_L;

model2.Q = sparse(eye(size_H));
model2.A = sparse([Aeq2; Af2]);
model2.rhs = [Beq2; theta0];
model2.sense = '=';
model2.quadcon(1).Qc = sparse(ones_H);
model2.quadcon(1).q = sparse(zeros(length(H_t)/4, 1));
model2.quadcon(1).rhs = umax^2;
model2.lb = -ub_L;
model2.ub = ub_L;

model3.Q = sparse(eye(size_H));
model3.A = sparse([Aeq3; Af3]);
model3.rhs = [Beq3; theta0];
model3.sense = '=';
model3.quadcon(1).Qc = sparse(ones_H);
model3.quadcon(1).q = sparse(zeros(length(H_t)/4, 1));
model3.quadcon(1).rhs = umax^2;
model3.lb = -ub_L;
model3.ub = ub_L;

model4.Q = sparse(eye(size_H));
model4.A = sparse([Aeq4; Af4]);
model4.rhs = [Beq4; theta0];
model4.sense = '=';
model4.quadcon(1).Qc = sparse(ones_H);
model4.quadcon(1).q = sparse(zeros(length(H_t)/4, 1));
model4.quadcon(1).rhs = umax^2;
model4.lb = -ub_L;
model4.ub = ub_L;



for k = 1:iterations
    model1.rhs = [Beq1; theta_k(:, k)];
    model2.rhs = [Beq2; theta_k(:, k)];
    model3.rhs = [Beq3; theta_k(:, k)];
    model4.rhs = [Beq4; theta_k(:, k)];
    
    result1 = gurobi(model1, params);
    result2 = gurobi(model2, params); 
    result3 = gurobi(model3, params);
    result4 = gurobi(model4, params);
    
    lambda1 = -result1.pi;
    lambda2 = -result2.pi;
    lambda3 = -result3.pi;
    lambda4 = -result4.pi;
    
    theta_k(:, k+1) = theta_k(:, k) + alpha_k * (lambda1(end-3:end) + lambda2(end-3:end) + lambda3(end-3:end) + lambda4(end-3:end));

end

theta_opt = x_opt1(end-3:end);
rel_err_con = zeros(1, iterations+1);
error = theta_k(:, :, 1) - theta_opt;
rel_err_con(1, :) = vecnorm(error) / norm(theta_opt);

figure(29)
semilogy(0:iterations, rel_err_con(1, :));
grid on;
title("Convergence of target position for quadratic constraint.");
xlabel('Iteration k');
ylabel('$$\frac{||\theta(k)-\theta^*||}{||\theta^*||}$$', 'Interpreter', 'Latex');

%% Question 2
options = optimoptions(@quadprog, 'Display', 'None');
theta0 = [0;0;0;0];
rho = [2 3 4 5 6 ];
iterations = 100;

theta_k = zeros(4, iterations+1, 4, length(rho));
for i = 1:4
    for j = 1:length(rho)
    theta_k(:, 1, i, j) = theta0;
    end
end

theta_avg = zeros(4, length(rho));

lambda_1 = zeros(4, length(rho));
lambda_2 = zeros(4, length(rho));
lambda_3 = zeros(4, length(rho));
lambda_4 = zeros(4, length(rho));

u_1 = zeros(size_u, 1);
u_2 = zeros(size_u, 1);
u_3 = zeros(size_u, 1);
u_4 = zeros(size_u, 1);

H1_ADMM = zeros(10, 10, length(rho));
H2_ADMM = zeros(10, 10, length(rho));
H3_ADMM = zeros(10, 10, length(rho));
H4_ADMM = zeros(10, 10, length(rho));

c1_ADMM = zeros(length(rho), 10);
c2_ADMM = zeros(length(rho), 10);
c3_ADMM = zeros(length(rho), 10);
c4_ADMM = zeros(length(rho), 10);

for j = 1:length(rho)
    H1_ADMM(:, :, j) = H1 + rho(j) * (BT1.' * BT1);
    H2_ADMM(:, :, j) = H2 + rho(j) * (BT2.' * BT2);
    H3_ADMM(:, :, j) = H3 + rho(j) * (BT3.' * BT3);
    H4_ADMM(:, :, j) = H4 + rho(j) * (BT4.' * BT4);

    c1_ADMM(j, :) = c1 + rho(j) * (x01.'*AT1.' - theta_avg(:, j).') * BT1 + lambda_1(:, j).'*BT1;
    c2_ADMM(j, :) = c2 + rho(j) * (x02.'*AT2.' - theta_avg(:, j).') * BT2 + lambda_2(:, j).'*BT2;
    c3_ADMM(j, :) = c3 + rho(j) * (x03.'*AT3.' - theta_avg(:, j).') * BT3 + lambda_3(:, j).'*BT3;
    c4_ADMM(j, :) = c4 + rho(j) * (x04.'*AT4.' - theta_avg(:, j).') * BT4 + lambda_4(:, j).'*BT4;
end

ulim = umax / Tfinal * ones(size_u, 1);
for k = 1:iterations
    for j= 1:length(rho)
        [u_opt1, fval,exitflag,output,lambda1] = quadprog(H1_ADMM(:, :, j), c1_ADMM(j, :), [], [], [], [], -ulim, ulim, [], options); 
        [u_opt2, fval,exitflag,output,lambda2] = quadprog(H2_ADMM(:, :, j), c2_ADMM(j, :), [], [], [], [], -ulim, ulim, [], options); 
        [u_opt3, fval,exitflag,output,lambda3] = quadprog(H3_ADMM(:, :, j), c3_ADMM(j, :), [], [], [], [], -ulim, ulim, [], options); 
        [u_opt4, fval,exitflag,output,lambda4] = quadprog(H4_ADMM(:, :, j), c4_ADMM(j, :), [], [], [], [], -ulim, ulim, [], options); 

        theta_k(:, k+1, 1, j) = AT1 * x01 + BT1 * u_opt1;
        theta_k(:, k+1, 2, j) = AT2 * x02 + BT2 * u_opt2;
        theta_k(:, k+1, 3, j) = AT3 * x03 + BT3 * u_opt3;
        theta_k(:, k+1, 4, j) = AT4 * x04 + BT4 * u_opt4;

        theta_avg(:, j) = mean(theta_k(:, k+1, :, j), 3);
        lambda_1(:, j) = lambda_1(:, j) + rho(j) * (theta_k(:, k+1, 1, j)-theta_avg(:, j));
        lambda_2(:, j) = lambda_2(:, j) + rho(j) * (theta_k(:, k+1, 2, j)-theta_avg(:, j));
        lambda_3(:, j) = lambda_3(:, j) + rho(j) * (theta_k(:, k+1, 3, j)-theta_avg(:, j));
        lambda_4(:, j) = lambda_4(:, j) + rho(j) * (theta_k(:, k+1, 4, j)-theta_avg(:, j));

        c1_ADMM(j, :) = c1 + rho(j) * (x01.'*AT1.' - theta_avg(:, j).') * BT1 + lambda_1(:, j).'*BT1;
        c2_ADMM(j, :) = c2 + rho(j) * (x02.'*AT2.' - theta_avg(:, j).') * BT2 + lambda_2(:, j).'*BT2;
        c3_ADMM(j, :) = c3 + rho(j) * (x03.'*AT3.' - theta_avg(:, j).') * BT3 + lambda_3(:, j).'*BT3;
        c4_ADMM(j, :) = c4 + rho(j) * (x04.'*AT4.' - theta_avg(:, j).') * BT4 + lambda_4(:, j).'*BT4;
    end
end

theta_opt = x_opt1(end-3:end);
rel_err = zeros(length(rho), iterations+1, 4);
for i = 1:length(rho)
    for j = 1:4
        error = theta_k(:, :, j, i) - theta_opt;
        rel_err(i, :, j) = vecnorm(error) / norm(theta_opt);
    end
end
for i = 1:4
    figure(29+i);
    for j = 1:length(rho)
        semilogy(0:iterations, rel_err(j, :, i)); hold on;
    end
    hold off;
    grid on;
    title(['Convergence of target position with ADMM for airplane ', num2str(i)]);
    xlabel('Iteration k');
    ylabel('$$\frac{||\theta(k)-\theta^*||}{||\theta^*||}$$', 'Interpreter', 'Latex');
    legend('\rho=2','\rho=3', '\rho=4', '\rho=5', '\rho=6'); 
end

figure(34);
semilogy(0:iterations, rel_err(3, :, 1)); 
grid on;
title('Convergence of target position with ADMM for airplane 1');
xlabel('Iteration k');
ylabel('$$\frac{||\theta(k)-\theta^*||}{||\theta^*||}$$', 'Interpreter', 'Latex');

xopt1 = findStates(AL1, BL1, AT1, BT1, x01, u_opt1);
xopt2 = findStates(AL2, BL2, AT2, BT2, x02, u_opt2);
xopt3 = findStates(AL3, BL3, AT3, BT3, x03, u_opt3);
xopt4 = findStates(AL4, BL4, AT4, BT4, x04, u_opt4);

n = length(A1);

for i = 1:n
    figure(34+i);
    plot(0:5, xopt1(i:n:end), '-o'); hold on;
    plot(0:5, xopt2(i:n:end), '-x'); 
    plot(0:5, xopt3(i:n:end), '-d'); 
    plot(0:5, xopt4(i:n:end), '-p'); hold off;
    legend('x_1', 'x_2', 'x_3', 'x_4');
    grid on;
    title(['State ', num2str(i)]);
    xlabel('Timestep t');
    ylabel('Position');
end



function [c, ceq] = nonlcon(x)

    c = [x(1:10).' * x(1:10) - 100^2;
         x(11:20).' * x(11:20) - 100^2;
         x(21:30).' * x(21:30) - 100^2;
         x(31:40).' * x(31:40) - 100^2];
    ceq = 0;
end
