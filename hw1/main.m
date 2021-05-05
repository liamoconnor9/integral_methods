%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calls the solve_fredholm function to answer the assignment questions
% and produce figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc;

% Plots error vs resolution for Trapezoid rule
if 0
    N_vec = [10:10:1000];
    k_vec = [1 4 5 15];
    labels = cell(length(k_vec), 1);
    count = 1;
    figure('units','normalized','outerposition',[0 0 1 1])
    for k = k_vec
        error_vec = 0 * N_vec;
        past = 0;
        for i = 1:length(N_vec)
            [x, t, e] = solve_fredholm(N_vec(i), k, 1);
            error_vec(i) = e;
            if (e < 1e-6 && ~past)
                past = 1;
                disp(strcat("For k = ", num2str(k), ", error<1e-6 at N = ", num2str(N_vec(i))));
            end
        end
        if k == k_vec(length(k_vec) - 1)
            loglog(N_vec, error_vec, '--', 'LineWidth', 5)
        else
            loglog(N_vec, error_vec, 'LineWidth', 5)
        end
        hold on
        labels{count} = strcat("k = ", num2str(k));
        count = count + 1;
        P = polyfit(log10(N_vec), log10(error_vec),1);
        disp(strcat("for k = ", num2str(k),  ", slope = ", num2str(P(1)), ", intercept = ", num2str(P(2))));
    end
    xlabel("$N$", 'interpreter', 'latex')
    ylabel("$|\mathbf{x} - \mathbf{x}_e|_{\infty}$", 'interpreter', 'latex')
    title("Trapezoidal Rule Error")
    legend(labels, 'Location', 'best')
    set(gca,'FontSize', 50)
    saveas(gcf, strcat('trap_error.png'))
    close
end

% Plots error vs resolution for Gauss-Legendre rule
if 0
    N_vec = [2:1:20];
    k_vec = [1 4 15];
    labels = cell(length(k_vec), 1);
    count = 1;
    figure('units','normalized','outerposition',[0 0 1 1])
    for k = k_vec
        error_vec = 0 * N_vec;
        past = 0;
        for i = 1:length(N_vec)
            [x, t, e] = solve_fredholm(N_vec(i), k, 0);
            error_vec(i) = e;
            if (e < 1e-6 && ~past)
                past = 1;
                disp(strcat("For k = ", num2str(k), ", error<1e-6 at N = ", num2str(N_vec(i))));
            end
        end
        if k == k_vec(length(k_vec) - 1)
            loglog(N_vec, error_vec, '--', 'LineWidth', 5)
        else
            loglog(N_vec, error_vec, 'LineWidth', 5)
        end
        hold on
        labels{count} = strcat("k = ", num2str(k));
        count = count + 1;
    end
    xlabel("$N$", 'interpreter', 'latex')
    ylabel("$|\mathbf{x} - \mathbf{x}_e|_{\infty}$", 'interpreter', 'latex')
    title("Gauss-Legendre Quadrature Error")
    legend(labels, 'Location', 'southwest' )
    set(gca,'FontSize', 50)
    saveas(gcf, strcat('gl_error.png'))
    close
end

N_vec = [1 2 3 4 5 6 7 8 9 10 11]
for N = N_vec
    [x, t, e] = solve_fredholm(N, 7, 1);
    disp(strcat("At N = ", num2str(N), ", max error = ", num2str(e)));
end