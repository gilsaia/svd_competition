function [u,s,v]=my_svd_1(A,r)
    [m,n]=size(A);
    [u,b,v]=bidiagonal_r(A,r,m,n);
    [u(:,1:n),s,v]=dk_svd(b,u(:,1:n),v,n);
    [s,v]=change_signval(s,v,r);
    s=[s;zeros(m-n)];
end