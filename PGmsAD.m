function [x_star, error] = PGmsAD(A, B, b)

%min_x max_{y, z} delta_{R+}(x) - x^T (A + B)^T y + b^T y - delta_{R+}(z)
% s.t. -x + (B - A)^T y + z = 0
%delta_{R+}为R+上的指示函数

%(y, z)为算法4.1里的y, (u, v)为算法4.1里的z

[n, ~] = size(A);
max_N = 10;
max_T = 1000;
alpha = 0.02;
alpha_x = alpha;
alpha_y = alpha;
alpha_z = alpha ;
 
x = zeros(n, 1);
x = [2.28571410e-01, 0.00000000e+00, 7.50000595e-09, 1.09523807e+01]';

y = zeros(n, 1);
y = (B - A)'\x;

z = zeros(n, 1);
u = zeros(n, 1);
u = y;
v = zeros(n, 1);


lambda = zeros(n, 1);



mu = 1;
tmp = max((B - A) \ (b - (A + B) * x), 0);
x_tmp = x - tmp;
error = norm(A * x_tmp + B * abs(x_tmp) - b, 2);

Q1 = @(x, y, z, lambda) b' * y - x' * (A + B)' * y - (mu/2) * (y - u)' * (y - u) + lambda' * (-x + (B - A)' *  y + z);

for iter = 1:max_T
    %内层循环(y, z)
    for k = 1:max_N
        grad_y = b - (A + B) * x + (B - A) * lambda - mu * (y - u);
        y = y + alpha_y * grad_y;
        grad_z = lambda - mu * (z - v);
        z = max(z + alpha_z * grad_z, 0);
        temp_y = Q1(x, y, z, lambda);
        disp(temp_y);
    end
    %temp_y = Q1(x, y, z, lambda)
    grad_x = -(A + B)' * y - lambda;
    u = (1 + alpha_x * mu) * u - alpha_x * mu * y;
    v = (1 + alpha_x * mu) * v - alpha_x * mu * z;
    lambda = lambda - alpha_x * (-x + (B - A)' * y + z);
    x = max(x - alpha_x * grad_x, 0);
    %temp_x = Q1(x, y, z, lambda)
    tmp = max((B - A) \ (b - (A + B) * x), 0);
    x_tmp = x - tmp;
    error(end + 1) = norm(A * x_tmp + B * abs(x_tmp) - b, 2);
end
x_star = x_tmp;



