clc; clear; close all;

clear all;

sizes = [1e3, 1e4, 1e5, 1e6]; % 测试不同规模的向量
num_trials = 1;             % 每个规模测试次数
n = 10000;
pl = @(x) (abs(x) + x)/2; 

% 预定义存储结果的变量
x_norms = zeros(num_trials, 1);
r_norms = zeros(num_trials, 1);
grad_norms = zeros(num_trials, 1);
iters = zeros(num_trials, 1);
time_newton = zeros(num_trials, 1);
time_prox = zeros(num_trials, 1);
x_init = zeros(n, 1);
r_norms_prox = zeros(num_trials, 1);
tol = 1e-23;
for k = 1:num_trials
    % 每次试验生成新的随机数据
    x = 10 * (rand(n, 1) - rand(n, 1));
    u = 1 * pl(x);
    A = null(u');
    A = A * A';
    A = A + 2 * eye(n);
    b = 5 * rand(n, 1);
    
    % 调用求解函数
    [x_star, r_star, iter, grad_norm, time1] = solve_AVE_GN(A, b);
    [x_star1, f_vals, time2] = solve_ave(A, b, x_init, tol);
    [x_star2, f_vals2, time3] = RIM(A, b, x_init, tol);
    
    % 记录当前试验结果
    x_norms(k) = norm(x_star, inf);
    r_norms(k) = norm(r_star, 2)^2;
    grad_norms(k) = grad_norm;
    iters(k) = iter;
    time_prox(k) = time2;
    time_newton(k) = time1;
    r_norms_prox(k) = f_vals(end);
end

% 计算平均值
avg_x_norm = mean(x_norms);
avg_r_norm = mean(r_norms);
avg_r_norm_prox = mean(r_norms_prox);
avg_grad_norm = mean(grad_norms);
avg_iter = mean(iters);
avg_time_newton = mean(time_newton);
avg_time_prox = mean(time_prox);
% 显示平均结果

disp(['Average time newton: ', num2str(avg_time_newton)]);
disp(['Average time prox: ', num2str(avg_time_prox)]);
disp(['Average r_norm newton: ', num2str(avg_r_norm)])
disp(['Average r_norm prox: ', num2str(avg_r_norm_prox)])
%{
disp(['Average x_norm: ', num2str(avg_x_norm)]);
disp(['Average r_norm: ', num2str(avg_r_norm)]);
disp(['Average grad_norm: ', num2str(avg_grad_norm)]);
disp(['Average iterations: ', num2str(avg_iter)]);

fprintf('%.7f\n', f_vals(end));
fprintf('%.7f\n', r_norms);
%}