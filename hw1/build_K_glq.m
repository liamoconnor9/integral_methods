%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Builds the kernel matrix and solves for the quadrature points 
%%% of the Gaussian Legendre Quadrature (glq) method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K, t_vec] = build_K_glq(N, k)

% Calculates Legendre coefficients at order N
coeff = LegendrePoly(N);

% These are the roots of the nominal Legendre Polynomial (on -1 < x < 1)
t_vec_nom = sort(roots(coeff));

% Our quadrature points take on the roots of the Nth order polynomial
% the scaling factor is pi / 4
t_vec = sort(pi * (t_vec_nom + 1) / 4.0);

% Calculates the coefficients for the Legendre Polynomial
pprime = polyder(coeff);

% Here we construct and populate a max (with congruent rows) with the appropriate 
% weights recall that the weights are obtained before the transformation and then
% muliplied by the scaling factor. The pi in the scaling factor cancels with the 1/pi in 
% K so we ignore them both
w_mat = zeros(N, N);
for col = 1:N
    w_mat(:, col) = (1/4)*(2/((1-(t_vec_nom(col))^2)*(polyval(pprime,t_vec_nom(col)))^2));
end

% Construct K and multiply element-wise with the weights matrix
K =  cos(k * (ones(N, 1) * t_vec' + t_vec * ones(1, N))) .* w_mat;

end