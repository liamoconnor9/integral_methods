%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Evaluates the the vector y at the given quadrature points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = y_func(t_vec, k)
    y = 0.75 * cos(k * t_vec) + sin(k * t_vec) * (1 - (-1)^k) / 4.0 / pi / k;
end