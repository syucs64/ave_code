function [x_star, f_vals, time] = solve_ave(A, b, x_init, tol)


max_iter = 10000;
x = x_init;
AT = A';
x_old = x;
grad_old = 2 * (AT * (A * x - b) + x);
L1 = 2 * (norm(A, 2)^2 + 1);
L2 = norm(A, 2);
f_vals = [sum((A * x - abs(x) - b).^2)];
r = A * x - b;
absx = abs(x);

tstart = tic;
for iter = 1:max_iter
    if iter == 1
        alpha = 1 / (L1 + 4 * L2);
    else
        delta_x = x - x_old;
        grad = 2 * (AT * r + x);
        delta_grad = grad - grad_old;
        m1 = sum(delta_x.^2);
        n1 = abs(delta_x' * delta_grad);
        alpha = m1 / n1 ;
        x_old = x;
        grad_old = grad; 
    end
    v1 = -2 * alpha * r;
    v2 = x - 2 * alpha * (AT * r + x - AT * absx);
    x_new = prox(v1, v2); 
    r_new = A * x_new - b;
    absx_new = abs(x_new);
    %%{
    ind = 1;
    while sum((r_new - absx_new).^2) > f_vals(end)- 1e-10 * sum((x - x_new).^2) && ind <= 10
        v1 = 0.5 * v1;
        v2 = 0.5 * (v2 + x);
        x_new = prox(v1, v2);
        r_new = A * x_new - b;
        absx_new = abs(x_new);
        ind = ind + 1;
    end
    %}
    x = x_new;
    r = r_new;
    absx = absx_new;
    f_vals(end + 1) = sum((r - absx).^2);
    if norm(f_vals(end) - f_vals(end - 1), 2) < tol
        break;
    end
end
time = toc(tstart);
x_star = x;

