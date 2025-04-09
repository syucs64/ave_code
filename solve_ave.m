function [x_star, f_vals] = solve_ave(A, b, x_init, tol)

n = size(A, 1);
max_iter = 10000;
x = x_init;
grad_h = @(x) A' * (A * x - b) + 2 * x;
f = @(x) norm(A * x - abs(x) - b, 2)^2;
x_old = x;
grad_old = grad_h(x_old);
L1 = 2 * (norm(A, 2)^2 + 1);
L2 = norm(A, 2);
f_vals = [f(x)];

for iter = 1:max_iter
    if iter == 1
        alpha = 1 / (L1 + 4 * L2);
    else
        delta_x = x - x_old;
        delta_grad = grad_h(x) - grad_old;
        m1 = norm(delta_x, 2)^2;
        n1 = abs(delta_x' * delta_grad);
        alpha = m1 / n1 ;
        x_old = x;
        grad_old = grad_h(x_old);
        %alpha = 1 / (L1 + 4 * L2);
    end
    v1 = -2 * alpha * (A * x - b);
    v2 = x - 2 * alpha * (A' * (A * x - b) + x - A' * abs(x));
    temp = prox(v1, v2); 
    const = norm(A * x - b, 2)^2 + norm(x, 2)^2;
    %g1 = -2 * (A * x - b)' * abs(x) + const;
    %g2 = (temp - x)' * (2  * (A' * (A * x - b) + x - A' * abs(x))) - 2 * (A * x - b)' * abs(temp) ...
        %+ 1 / (2 * alpha) * (temp - x)' * (temp - x) + const;
    x = temp;
    f_vals(end + 1) = f(x);
    if norm(x - x_old, 2) < tol
        break;
    end
end

x_star = x;

