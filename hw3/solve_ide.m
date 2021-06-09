%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solves the integro-differential equation:
% $\frac{d^2h}{dx^2} + h + \frac{i}{2} \int_{-\pi/2}^{\pi/2} H^1_0 \big( |x - y| \big)...
% h(y)dy &= e^{iax}$
% via finite difference, trapezoidal rule, product integration, linear
% interpolation, and finally some direct solve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;
N = 128;
a = 0.1;
% filename = 'ide21';

pts = 7;
error = zeros(pts - 1);
res = zeros(pts - 1);
for i = 1:pts
    n = 2 * 2^i;
    [h, x] = solve_ide(a, n);
    h(1) = 0.0;
    h(N + 1) = 0.0;
 
    if (i == 1)
        ho = max(abs(h));
    else
        error(i - 1) = abs(max(abs(h)) - ho);
        res(i - 1) = n;
        disp(abs(max(abs(h)) - ho));
        ho = max(abs(h));
    end
end

figure('units','normalized','outerposition',[0 0 1 1])
plot(res, error, 'LineWidth', 5);
% hold on
% plot(x, imag(h), '--', 'LineWidth', 5);
xlabel("$N$", 'interpreter', 'latex')
ylabel("$|h_{N} - h_{N/2}|_{\infty}$", 'interpreter', 'latex')
title(strcat("Error Difference: a = ", num2str(a)))
% legend('Real', 'Imaginary', 'Location', 'best')
set(gca,'FontSize', 50)
saveas(gcf, strcat('errors1', '.png'))
close


% figure('units','normalized','outerposition',[0 0 1 1])
% plot(x, real(h), 'LineWidth', 5);
% hold on
% plot(x, imag(h), '--', 'LineWidth', 5);
% xlabel("$x$", 'interpreter', 'latex')
% ylabel("$h(x)$", 'interpreter', 'latex')
% title(strcat("Equation 1 Solution: a = ", num2str(a)))
% legend('Real', 'Imaginary', 'Location', 'best')
% set(gca,'FontSize', 50)
% saveas(gcf, strcat(filename, '.png'))
% close



function [h, x_vec] = solve_ide(a, N)

delx = pi / N;
x_vec = linspace(-pi/2, pi/2, N + 1);

% 2nd derivative finite difference matrix 
D = (-2*diag(ones(1,N+1)) + diag(ones(1,N),1) + diag(ones(1, N), -1)) / delx / delx;
I = eye(N+1);

% trapezoidal rule weights
w = delx * ones(N+1); 
w(1) = w(1) / 2;
w(N+1) = w(N+1) / 2;

% conventional trapezoidal kernels for the nonsingular terms
% Bessel function
K1 = zeros(N+1);
for i = 1:N+1
    for j = 1:N+1
        K1(i, j) = w(j) * besselj(0, abs(x_vec(i) - x_vec(j)));
    end
end

% Polynomial (infinite series) term associated with the Bessel function of
% the second kind Y0(z)
K2 = zeros(N+1);
for i = 1:N+1
    for j = 1:N+1
        dx = abs(x_vec(i) - x_vec(j));
        for k = 0:20
                % Explicit equations for the weights given by Viktorovitch et al
                % https://www.witpress.com/Secure/elibrary/papers/BE95/BE95025FU.pdf
                % in Appendix A
                % They report a tolerance of 1.4e-8 in 0<x<3 using on 7
                % terms. The magnitude of the weights drops off rapidly so,
                % at some point, it's not product to keep adding terms!
                weight = -(-1)^k * psi(k + 1) / pi / (4^k) / factorial(k) / factorial(k);
                K2(i, j) = K2(i, j) + w(j) * weight * dx ^ (2 * k);
        end
    end
end


% Here we must approximate a more complicated integral 
% to obtain the weights which correspond to the linear problem
A = zeros(N+1);
for i = 1:N+1
    % The rows of A correspond to the quadrature points
    for j = 2:N
        f_neg = @(x) (x - x_vec(j-1)) .* besselj(0, abs(x_vec(i) - x)) .* log((abs(x_vec(i) - x)) / 2);
        f_pos= @(x) (x_vec(j+1) - x) .* besselj(0, abs(x_vec(i) - x)) .* log((abs(x_vec(i) - x)) / 2);
        A(i, j) = integral(f_neg, x_vec(j-1), x_vec(j));
        A(i, j) = A(i, j) + integral(f_pos, x_vec(j), x_vec(j+1));
        A(i, j) = A(i, j) / delx / pi;
    end
end

M = vertcat([D + I - (K2 + A) , -K1 / 2], [K1 / 2 , D + I - (K2 + A) / 2]);
b = vertcat(cos(a * x_vec'), sin(a * x_vec'));
x = linsolve(M, b);
h = x(1:N + 1) + 1i * x(N + 2:2*N + 2);

end