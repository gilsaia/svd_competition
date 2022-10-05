function [u,s,v]=my_svd_1(A,r)
    [m,n]=size(A);
    [u,b,v]=bidiagonal_r(A,r,m,n);
    % [ux,s,v]=jacobi(b,v,n,r);
    % s=[s;zeros(m-n)];
    % u(:,1:n)=matmul(u(:,1:n),ux);
    [u(:,1:n),s,v]=dk_svd(b,u(:,1:n),v,n);
    s=[s;zeros(m-n)];
end