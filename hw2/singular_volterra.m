%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solves the Singular Volterra Integral Equation:
% $x(s) = (1 + s)^{-1/2} + \frac{\pi}{8} - \frac{1}{4}\arcsin \Big( \frac{1 - s}{1 + s} \Big) 
%  - \frac{1}{4} \int_{0}^{s} \frac{x(t)}{(s - t)^{1/2}} dt$ on $0 < s < 1$
%   via Trapezoidal product integration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters
N = 16;

% Allocation
s = linspace(0, 1, N + 1)';
se = linspace(0, 1, 1000)';
b = 1 ./ sqrt(1 + s) + pi / 8 - asin((1 - s) ./ (1 + s)) / 4;
W = zeros(N + 1, N + 1);
xe = 1 ./ sqrt(s + 1);
xe_plot = 1 ./ sqrt(se + 1);

% We let the first row have zero weights because this is associated with  
% s = 0. The contribution of the integral term is zero here so the weight
% shall be as well.
for row = 2:N + 1 %i
    for col = 1:row %j
        if col == 1
            W(row, col) = (row - 2)^(3/2) - sqrt(row - 1) * (2 * (row - 1) - 3)/2;
        elseif col == row
            W(row, col) = 1;
        else
            W(row, col) = (row - col + 1)^(3/2) - 2 * (row - col)^(3/2) + (row - col - 1)^(3/2);
        end
    end
end

x = linsolve(eye(N + 1) + W / sqrt(N)/3,  b);

figure('units','normalized','outerposition',[0 0 1 1])
plot(s, x, 'LineWidth', 5);
hold on
plot(se, xe_plot, '--', 'LineWidth', 5);
xlabel("$s$", 'interpreter', 'latex')
ylabel("$x$", 'interpreter', 'latex')
title(strcat("Volterra Equation (1) Solutions: N = ", num2str(N)))
legend('Numerical', 'Analytical', 'Location', 'best')
set(gca,'FontSize', 50)
saveas(gcf, strcat('volt', num2str(log2(N)), '.png'))
close