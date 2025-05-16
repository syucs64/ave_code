clc; clear;
n = 100;
I_n = -eye(n);
    
% 构建右上角的 n×2 零矩阵
zero_block = zeros(n, 2);

% 构建上半部分 [ -I_n | zeros(n,2) ]
upper_part = [I_n, zero_block];

% 构建 e_n^T (全1的行向量)
e_n = ones(1, n);

% 构建左下角的 2×n 矩阵 [ e_n^T; -e_n^T ]
lower_left = [e_n; -e_n];

% 构建右下角的 2×2 矩阵 [ -1 0; 0 -1 ]
lower_right = [-1 0; 0 -1];

% 构建下半部分 [ e_n^T -1 0; -e_n^T 0 -1 ]
lower_part = [lower_left, lower_right];

% 组合成完整的矩阵
M = [upper_part; lower_part];
I = eye(n + 2);
A = (M - I) \ (M + I);

[U, S, V] = svd(A);
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