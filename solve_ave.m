function [x_star, f_vals] = solve_ave(A, b, x_init, tol)

n = size(A, 1);
max_iter = 1000;
x = x_init;
grad_h = @(x) A' * (A * x - b) - 2 * x;
f = @(x) norm(A * x - abs(x) - b, 2);
x_old = x;
grad_old = grad_h(x_old);
L1 = 2 * (norm(A, 2)^2 + 1);
L2 = norm(A, 2);
f_vals = [];

for iter = 1:max_iter
    if iter == 1
        alpha = 1 / (L1 + 4 * L2);
    else
        delta_x = x - x_old;
        delta_grad = grad_h(x) - grad_old;
        m = norm(delta_x, 2)^2;
        n = abs(delta_x' * delta_grad);
        alpha = m / n ;
        x_old = x;
        grad_old = grad_h(x_old);
    end
    temp = zeros(n, 1);
    v1 = -2 * alpha * A * x - b;
    v2 = x - 2 * alpha * (A' * (A * x - b) + x - A' * abs(x));
    for i = 1:n
        a = v1(i);
        u = v2(i);
        temp(i) = prox(a, u);
    end
    x = temp;
    f_vals(end + 1) = f(x);
    if norm(x - x_old, 2) < tol
        break;
    end
end

x_star = x;

