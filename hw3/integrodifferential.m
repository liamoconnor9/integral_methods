       % To solve: the following Integro Differential Equation
% h''(x) + h(x) + (i/2) * Int[-Pi/2 to Pi/2] H_0^1(|x-y|)h(y)dy = e^{iax}
                    %h(Pi/2) = h(-Pi/2) = 0

                    % h(x) = hr + i hi(x)
                    
function [u, ur, ui] = integrodifferential(a, N)

I = eye(N+1); %Identity Matrix

h = pi/N; %Step Size

%Central Difference Second Derivative Operator
L = diag(-2*(1/h^2)*ones(1,N+1)) + diag((1/h^2)*ones(1,N),1) + diag((1/h^2)*ones(1,N),-1);

%Quadrature points
s = -pi/2:h:pi/2; 

%Operators%

%Repeated Trapezoidal Rule Weights
w = h.*ones(N+1); 
w(1) = w(1)/2;
w(N+1) = w(N+1)/2;

%Bessel first kind operator
K1 = zeros(N+1);
for i = 1:N+1
    for j = 1:N+1
        K1(i,j) = w(j)*besselj(0,abs(s(i)-s(j)));
    end
end

%Regular part of Webber operator
K2 = zeros(N+1);
for i = 1:N+1
    for j = 1:N+1
        dij = abs(s(i)-s(j))/3;
        K2(i,j) = w(j)*(0.36746691 + 0.60559366*dij^2 - 0.74350384*dij^4 ... 
                  + 0.25300117*dij^6 - 0.04261214*dij^8 + 0.00427916*dij^10 ...
                  - 0.00024846*dij^12);
    end
end

%Irregular part of Webber operator
K3 = zeros(N+1);
for i = 1:N+1
    for j = 2:N
        f1 = @(x) (2/pi).*(x-s(j-1)).*besselj(0, abs(s(i)-x)).*log((abs(s(i)-x))/2);
        f2 = @(x) (2/pi).*(s(j+1)-x).*besselj(0, abs(s(i)-x)).*log((abs(s(i)-x))/2);
        K3(i,j) = (1/h)* (integral(f1, s(j-1), s(j)) + integral(f2, s(j), s(j+1)));
    end
end

%Concatenated Matrix
A = [L+I-(K2+K3)/2 , -K1/2 ; K1/2 , L+I-(K2+K3)/2];

%Concatenated RHS Vector
fc = cos(a*s);
fs = sin(a*s);
f = [fc, fs];

%Solving for u
u = A\f';

%Extracting real and complex parts from the solution list
ur = u(1:N+1);
ui = u(N+2:2*N+2);

%The final complex-valued solution 
u = ur + 1i*ui;

end