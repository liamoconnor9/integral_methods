%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Builds the kernel matrix for repeated Trapezoid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K, t_vec] = build_K_trap(N, k)

% Uniform grid size
dt = pi / 2.0 / N;

% Quadrature points vector
t_vec = linspace(0, pi / 2.0, N)';

% K is almost symetric. We halve the endpoints but all other weights are one 
K = (dt * cos(k * (ones(N, 1) * t_vec' + t_vec * ones(1, N)))) / pi;
K(:,1) = K(:,1)/2;
K(:,N) = K(:,N)/2;

end