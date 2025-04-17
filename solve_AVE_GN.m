function [x_star, r_star, iter, grad_norm, time] = solve_AVE_GN(A, b, x_init,max_iter, tol, delta, s, mu)


if nargin < 4, max_iter = 1000; end
if nargin < 5, tol = 1e-12; end
if nargin < 6, delta = 1e-4; end
if nargin < 7, s = 1; end
if nargin < 8, mu = 0.5; end

n = size(A, 1);
x = x_init; % 初始点x0=0
iter = 0;
grad_norm = inf;
tstart = tic;
while iter < max_iter
    % 计算当前Q矩阵
    D = diag(sign(x));
    Q = A - D;
    
    % 计算梯度和Hessian
    residual = Q * x - b;
    grad = 2 * Q' * residual;
    H = 2 * (Q' * Q);
    
    % 检查收敛条件
    grad_norm = norm(grad, inf);
    if grad_norm < tol
        break;
    end
    
    % 计算修正的Hessian和牛顿方向
    H_mod = H + delta * eye(n);
    d = -H_mod \ grad;
    
    % Armijo线搜索
    alpha = s;
    armijo_iter = 0;
    max_armijo = 8;
    while armijo_iter < max_armijo
        x_new = x + alpha * d;
        D_new = diag(sign(x_new));
        Q_new = A - D_new;
        residual_new = Q_new * x_new - b;
        f_new = norm(residual_new)^2;
        f_current = norm(residual)^2;
        
        % Armijo条件
        if f_current - f_new >= -mu * alpha * grad' * d
            break;
        else
            alpha = alpha * 0.5; % 减小步长
            armijo_iter = armijo_iter + 1;
        end
    end
    
    % 更新x
    x = x + alpha * d;
    iter = iter + 1;
end

x_star = x;
r_star = A * x_star - abs(x_star) - b; % 根据式(7)计算r*
time = toc(tstart);
end