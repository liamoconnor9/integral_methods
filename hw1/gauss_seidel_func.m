%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solves the system Ax = b to some tolerance using Gauss Seidel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = gauss_seidel_func(A, x, b)

tol = 1e-5;

% Fixed parameters
N = length(A);
D = diag(diag(A));
L = tril(D - A);
DmL_inv = (D - L)^(-1);
G =  DmL_inv * transpose(L);

% residuals / error
r = norm(A*x - b, inf);
error_vec = [r];

% iterative operation
while r >= tol
    x = G * x + DmL_inv * b;
    r = norm(A*x - b, inf);
    error_vec = [error_vec, r];
    if (r > error_vec(end) || length(error_vec) > 10000)
        r = 0;
        disp("BAD BAD")
        throw
    end
%     disp(length(error_vec))
%     disp(r)
end