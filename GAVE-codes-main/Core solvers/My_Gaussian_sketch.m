function [SA, SB, Sb] = My_Gaussian_sketch(A, B, b, ell)
    % Gaussian Sketch
    % Input:
    %   A: m x n matrix
    %   B: m x n matrix
    %   b: m x 1 vector
    %   ell: target dimension (ell < m)
    % Output:
    %   SA: ell x n compressed matrix
    %   SB: ell x n compressed matrix
    %   Sb: ell x 1 compressed vector
    % Coded by Jiaxin Xie, Beihang University, xiejx@buaa.edu.cn

    % Get the size of matrix A
    [m, n] = size(A);  

    % Check if ell is less than m
    if ell >= m
        error('Target dimension ell must be less than m.');
    end

    % Step 1: Generate a Gaussian random matrix S
    % S is an ell x m matrix with entries drawn from N(0, 1)
    S = randn(ell, m) / sqrt(ell);  

    % Step 2: Compute the Gaussian sketch for A, B, and b
    % Compute the sketched matrices and vector
    SA = S * A;  
    SB = S * B;  
    Sb = S * b;  
end
