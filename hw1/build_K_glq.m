function x = build_K_glq(N, k)

y = @(t) ((3/4)*cos(k*t)) + ((sin(k*t)*(1-(-1)^k))/(4*pi*k));
xexact = @(t) cos(k*t);

lcoef = LegendrePoly(N);
ldcoef = polyder(lcoef);

I = eye(N);


if methodid == 0
    dx = pi/(2*(N-1));
    svec = zeros(1, N);
    for i = 1:N
        svec(i) = (i-1)*dx;
    end
    K = zeros(N,N);
    for i=1:N
        for j = 1:N
            K(i,j) = (dx * cos(k*(svec(i)+svec(j))))/pi;
        end
    end
    K(:,1) = K(:,1)/2;
    K(:,N) = K(:,N)/2;
    
elseif methodid == 1
        svec1 = roots(lcoef);
        svec1 = sort(svec1);
        svec1 = reshape(svec1, 1, N);
        svec = (pi/4)*(svec1 + 1);
        w = zeros(N, 1);
        for i=1:N
            w(i) = (pi/4)*(2/((1-(svec1(i))^2)*(polyval(ldcoef,svec1(i)))^2));
        end
        K = zeros(N,N);
        for i=1:N
            for j=1:N
                K(i,j) = (w(j)*cos(k*(svec(i)+svec(j))))/pi;
            end
        end
end

A = I-K;

f = y(svec);
f = reshape(f, N, 1);

x = ones(N,1);
x = gs(A, f, 1e-6, x);

xex = xexact(svec);
xex = reshape(xex, N, 1);

error = norm(x-xex, inf);