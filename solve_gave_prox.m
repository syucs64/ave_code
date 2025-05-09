function [x_star, f_vals, time] = solve_gave_prox(A, B, b, x_init, tol)

% min ||Ax - B|x| - b||^2
max_iter = 100000;
x = x_init;
AT = A';
BT = B';
L1 = 2 * (norm(A, 2)^2);
L2 = norm(B, 2)^2;
L3 = 2 * norm(BT * A, 2);
f_vals = [sum((A * x - B * abs(x) - b).^2)];
absx = abs(x);
r = A * x -  B * absx - b;

alpha = 1 / (L1 + 2 * L2 + 2 * L3);
tstart = tic;
for iter = 1:max_iter
    v1 = -2 * alpha * BT * r;
    v2 = x - 2 * alpha * AT * r;
    x_new = prox(v1, v2); 
    absx_new = abs(x_new);
    r_new = A * x_new - B * absx_new - b;
    %%{
    ind = 1;
    while sum(r_new.^2) > f_vals(end)- 1e-10 * sum((x - x_new).^2) && ind <= 10
        v1 = 0.5 * v1;
        v2 = 0.5 * (v2 + x);
        x_new = prox(v1, v2);
        absx_new = abs(x_new);
        r_new = A * x_new - B * absx_new - b;
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

