function [x_star, f_vals] = PGmsAD(A, B, b)

[n, ~] = size(A);
max_N = 5;
max_T = 1000;
alpha_x = 0.02;
alpha_y = 0.02;
alpha_z = 0.02;
 
x = zeros(n, 1);
y = zeros(n, 1);
z = zeros(n, 1);
u = zeros(n, 1);
v = zeros(n, 1);
lambda = zeros(n, 1);

mu = 1;
tmp = max((B - A) \ (b - (A + B) * x), 0);
x_tmp = x - tmp;
f_vals = norm(A * x_tmp + B * abs(x_tmp) - b, 2);

Q1 = @(x, y, z, lambda) b'*y - x'*(A + B)'*y - (mu/2)*(y - u)'*(y - u) + lambda'*(-x + (B - A)'*y + z);
Q2 = @(x, y, z, lambda) -(mu/2)*(z - v)'*(z - v) + lambda'*(-x + (B - A)'*y + z);
Q3 = @(x, y, z, lambda) lambda'*(-x + (B - A)'*y + z) - x'*(A + B)'*y;
for iter = 1:max_T
    for k = 1:max_N
        ay = alpha_y;
        grad_y = b - (A + B) * x + (B - A) * lambda - mu * (y - u);
        y_tmp = y + ay * grad_y;
        while Q1(x, y_tmp, z, lambda) < Q1(x, y, z, lambda)
            ay = 0.5 * ay;
            y_tmp = y + ay * grad_y;
        end
        y = y_tmp;
        az = alpha_z;
        grad_z = lambda - mu * (z - v);
        z_tmp = z + az * grad_z;
        z_tmp = max(z_tmp, 0);
        while Q2(x, y, z_tmp, lambda) < Q2(x, y, z, lambda)
            az = 0.5 * az;
            z_tmp = z + az * grad_z;
            z_tmp = max(z_tmp, 0);
        end
        z = z_tmp;
    end
    ax = alpha_x;
    grad_x = -(A + B)' * y - lambda;
    x_new = x - ax * grad_x;
    x_new = max(x_new, 0);
    while Q3(x_new, y, z, lambda) > Q3(x, y, z, lambda)
        ax = 0.5 * ax;
        x_new = x - ax * grad_x;
        x_new = max(x_new, 0);
    end
    u = (1 + ax * mu) * u - ax * mu * y;
    v = (1 + ax * mu) * v - ax * mu * z;
    lambda = lambda - ax * (-x + (B - A)' * y + z);
    x = x_new;
    
    tmp = max((B - A) \ (b - (A + B) * x), 0);
    x_tmp = x - tmp;
    f_vals(end + 1) = norm(A * x_tmp + B * abs(x_tmp) - b, 2);
end
x_star = x_tmp;


