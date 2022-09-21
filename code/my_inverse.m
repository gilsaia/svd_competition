function [inv_AA] = my_inverse(A)
    [m,n]=size(A);
    I=eye(m,m);
    inv_AA=inv(I+A*A');
end