function [x_star, f_vals, time] = RIM(A, b, x_init, tol)

x = x_init;
AT = A';
A2 = norm(A, 2)^2;
maxiter = 10000;
A_pinv = pinv(A);
f_vals = [norm(A * x - abs(x) - b, 2)^2];
tstart = tic;
for i = 1:maxiter
    alpha = 0.5;
    x = x - alpha / A2 * AT * (A * x - abs(x) - b);
    %x = x - alpha * A_pinv * (A * x - abs(x) - b);
    %D = diag(sign(x));
    %x = (A - D)\ b;
    %x = A_pinv * (abs(x) + b);
    f_vals(end + 1) = sum((A * x - abs(x) - b).^2);
    if norm(f_vals(end) - f_vals(end - 1), 2) < tol
        break;
    end
end
time = toc(tstart);
x_star = x;

