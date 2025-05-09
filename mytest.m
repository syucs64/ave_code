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