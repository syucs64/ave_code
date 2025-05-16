%% run_minmax_AB.m  —— 直接运行的脚本
%   min_{x>=0}  g(x),   g(x)=max_y (b-Px)'y  s.t.  M' y <= x
%   A, B, b 见题；M=B-A, P=A+B
clear; clc;

%% ----------------- 数据 ------------------------------------------------
A = [  5     4   -3       0 ;
       4    11   -1     0.5;
      -3    -1    6     0.5;
       0   0.5  0.5   0.25];

B = [  2.5   1   -2       0 ;
        1   0.5  0.5  -0.25;
       -2   0.5  0.5  -0.25;
        0  -0.25 -0.25 -0.125];

b  = ones(4,1);
P  = A + B;
M  = B - A;   Mt = M.';  %#ok<NASGU>

%% ----------------- 外层优化设置 ---------------------------------------
x0  = 10*ones(4,1);      % **较大的正初值** ⇒ 内层先有界
lb  = zeros(4,1);        % x >= 0
ub  = 100*ones(4,1);     % 给一个比较宽松的上界

opts = optimoptions('fmincon',...
        'Algorithm','interior-point',...
        'Display','iter',...
        'FiniteDifferenceType','central',...
        'MaxFunctionEvaluations',2e4);

[xOpt,fOpt] = fmincon(@(x) outerObjective(x,b,P,M), ...
                      x0,[],[],[],[],lb,ub,[],opts);

%% 取对应 y*, z* 便于展示
[~,yOpt,zOpt] = innerMax(xOpt,b,P,M);

fprintf('\n===========  计算结果  ===========\n');
fprintf('g(x*) = %.12g\n', fOpt);
fprintf('x* = [%g  %g  %g  %g]^T\n', xOpt);
fprintf('y* = [%g  %g  %g  %g]^T\n', yOpt);
fprintf('z* = [%g  %g  %g  %g]^T\n', zOpt);

%% #######################################################################
%                   子  函  数
% #######################################################################
function [g,yStar,zStar] = innerMax(x,b,P,M)
    %  内层 max_y LP；若无界/不可行 → 返回大罚值
    Mt     = M.';
    q      = b - P*x;              % 目标 (b-Px)'y
    c      = -[q ; -q];            % linprog 最小化 —— 取负号
    Aineq  = [Mt, -Mt];            % y = y+ - y-, 约束 M'y <= x
    bineq  =  x;
    lb     = zeros(8,1);           % y+, y- >= 0
    BIG    = 1e20;                 % 罚值

    try
        optsLP = optimoptions('linprog','Display','none');
        [var,fval,exitflag] = linprog(c,Aineq,bineq,[],[],lb,[],optsLP);
    catch
        exitflag = -99;  % linprog 抛异常
    end

    if exitflag == 1                % 正常最优
        yStar = var(1:4) - var(5:8);
        zStar = x - Mt*yStar;
        g     = -fval;              % 恢复最大值
    else                            % 无界 / 不可行 / 异常
        yStar = [];  zStar = [];
        g     = BIG;
    end
end

function g = outerObjective(x,b,P,M)
    %  fmincon 只关心标量值
    g = innerMax(x,b,P,M);
end
