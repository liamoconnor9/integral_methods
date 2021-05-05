%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solves the system Ax = b to some tolerance using SOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, error_vec] = sor_func(A, x, b, w, tol)

% Fixed parameters
N = length(A);
I = eye(N);
D = diag(diag(A));
L = tril(D - A);
DmwL_inv = (D - w * L)^(-1);
Gw =  DmwL_inv * (w * transpose(L) + (1 - w) * D);
% disp(max(abs(eig(Gw))))

% residuals / error
r = norm(A*x - b, inf);
error_vec = [r];

% iterative operation
while r >= tol
    x = Gw * x + w * DmwL_inv * b;
    r = norm(A*x - b, inf);
    error_vec = [error_vec, r];
    if (r > error_vec(end) || length(error_vec) > 10000)
        r = 0;
    end
    %     disp(length(error_vec))
    %     disp(r)
end