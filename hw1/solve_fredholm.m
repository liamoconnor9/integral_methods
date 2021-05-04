solve_fredholm(100, 1, 1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Solves the Fredholm Equation of the second kind
%%% useTrap: 1 for Trapezoid method; 0 for Gaussian-Legendre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, t_vec, error] = solve_fredholm(N, k, useTrap)

    if (useTrap)
        [K, t_vec] = build_K_trap(N, k);
    else
        [K, t_vec] = build_K_glq(N, k);
    end
    
    y = y_func(t_vec, k);
    x = linsolve(eye(N) - K, y);
    
    plot(t_vec, x)
    
    error = 0;


end