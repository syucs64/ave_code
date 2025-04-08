clc; clear; close all;
n = 1000;
pl = @(x) (abs(x) + x)/2;  
num_trials = 50; % 设置随机试验次数

% 预定义存储结果的变量
x_norms = zeros(num_trials, 1);
r_norms = zeros(num_trials, 1);
grad_norms = zeros(num_trials, 1);
iters = zeros(num_trials, 1);

for k = 1:num_trials
    % 每次试验生成新的随机数据
    x = 10 * (rand(n, 1) - rand(n, 1));
    u = 1 * pl(x);
    A = null(u');
    A = A * A';
    A = A + eye(n);
    b = 5 * rand(n, 1);
    
    % 调用求解函数
    [x_star, r_star, iter, grad_norm] = solve_AVE_GN(A, b);
    
    % 记录当前试验结果
    x_norms(k) = norm(x_star, inf);
    r_norms(k) = norm(r_star, inf);
    grad_norms(k) = grad_norm;
    iters(k) = iter;
end

% 计算平均值
avg_x_norm = mean(x_norms);
avg_r_norm = mean(r_norms);
avg_grad_norm = mean(grad_norms);
avg_iter = mean(iters);

% 显示平均结果
disp(['Average x_norm: ', num2str(avg_x_norm)]);
disp(['Average r_norm: ', num2str(avg_r_norm)]);
disp(['Average grad_norm: ', num2str(avg_grad_norm)]);
disp(['Average iterations: ', num2str(avg_iter)]);