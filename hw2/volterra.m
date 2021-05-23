 % To solve: the following Integral Equation
% x(s)=(1+s)^(-1/2) + pi/8 - arcsin((1-s)/(1+s))/4 - (1/4)*Int[0 to s](s-t)^(-1/2)x(t)dt

%function [x, error] = volterra(N)

N = 30;

y = @(s) (1/sqrt(1+s)) + (pi/8) - (asin((1-s)/(1+s))/4); %The RHS function
xexact = @(s) 1/sqrt(1+s); %The exact solution

I = eye(N+1); %Identity Matrix

h = 1/N;

stencil = linspace(0,1,N+1);

%Defining the quadrature matrix
K = zeros(N+1);

for i = 2:N+1
    for j = 1:i
        if j == 1
            K(i, j) = (sqrt((i-2)*h)*(i-2))/3 - (sqrt((i-1)*h)*(2*i-5))/6;
        elseif j == i
            K(i, j) = sqrt(h)/3;
        else
            K(i, j) = (sqrt(h)/3)*((i-j-1)^(3/2)-2*(i-j)^(3/2)+(i-j+1)^(3/2));
        end
        
    end
end

%Solving the Linear System X = f + KX

A = I+K;

f = zeros(N+1, 1);

for i=1:N+1
    f(i) = y(stencil(i));
end

%x = ones(N+1,1);
%x = gs(A, f, 1e-6, x); %Gauss-Seidel

x = A\f;

xex = zeros(N+1, 1);

for i=1:N+1
    xex(i) = xexact(stencil(i));
end
plot(x)
error = norm(x-xex, inf); %Infinity Norm Error