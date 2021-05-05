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
%     [x, er] = sor_func(eye(N) - K, t_vec*0, y, 1.1, 1e-9);
    
%     error = norm(, inf);
er_vec = cos(k*t_vec) - x;
emin = min(er_vec);
emax = max(er_vec);
if (abs(emin) > abs(emax))
    error = emin;
else 
    error = emax;
end