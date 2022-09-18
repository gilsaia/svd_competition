function [C] = matmul(A, B)
% MATMUL   Matrix multiplication (product) implementation via definition.
%
%                C(i, j) = Sum[A(i, k) * B(k, j)]
%
%   C = MATMUL(A, B) returns the matrix product of A and B.
%   Where, the number of columns of A must be equal to the number of 
%   rows of B. When A is a row vector and B is a column vector, C will be a number.
%   
%   % Example:
%   %   matrix multiplication of two matrices.
%
%   rng default;
%   M = 1024; N = 512;
%   A = randn(M, N);
%   B = randn(N, M);
%   %  built in
%   C1 = A * B;
%   %  matmul
%   C2 = matmul(A, B);
%
%   See also EBEMUL.
%


narginchk(2,2)

[Ma, Na] = size(A);
[Mb, Nb] = size(B);

if Na ~= Mb
    error("The number of columns of A must be equal to the number of rows of B!");
end

C = zeros(Ma, Nb);
for i=1:Ma
    for j=1:Nb
        c = 0;
        for k=1:Na
             c = c + A(i, k) * B(k, j);
        end
        C(i, j) = c;
    end
end

end

