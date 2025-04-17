function [x_star, f_vals, time] = solve_ave_prox3(A, b, x_init, tol)

% min ||x - A^(-1)|x| - A^(-1)b||^2
max_iter = 10000;
x = x_init;
B = inv(A);
f = B * b;
L1 = norm(B, 2);
L2 = L1 ^ 2;
alpha = 1 / (4 * L1 + 2 * L2 + 1);
BT = B';
BTB = BT * B;
absx = abs(x);
r = x - B * absx - f;
f_vals = [sum(r.^2)];
tstart = tic;
for iter = 1:max_iter
    v1 = 2 * alpha * (BTB * absx - BT * (x - f));
    v2 = x - 2 * alpha * r;
    x_new = prox(v1, v2); 
    absx_new = abs(x_new);
    r_new = x_new - B * absx_new - f;
    %%{
    ind = 1;
    while sum(r_new.^2) > f_vals(end)- 1e-6 * sum((x - x_new).^2) && ind <= 10
        v1 = 0.5 * v1;
        v2 = 0.5 * (v2 + x);
        x_new = prox(v1, v2);
        absx_new = abs(x_new);
        r_new = x - B * absx_new - f;
        ind = ind + 1;
    end
    %}
    x = x_new;
    r = r_new;
    absx = absx_new;
    f_vals(end + 1) = sum(r.^2);
    if norm(f_vals(end) - f_vals(end - 1), 2) < tol
        break;
    end
end
time = toc(tstart);
x_star = x;

