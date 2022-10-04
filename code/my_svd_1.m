function [u,s,v]=my_svd_1(A,r)
    [m,n]=size(A);
    [u,b,v]=bidiagonal_r(A,r,m,n);
    [ux,s,v]=jacobi(b,v,n);
    s=[s;zeros(m-n)];
    uo=[ux zeros(m-n);zeros(m-n) eye(m-n)];
    u=matmul(u,uo);
    % [u(:,1:n),s,v]=dk_svd(b,u(:,1:n),v,n);
    % s=[s;zeros(m-n)];
end