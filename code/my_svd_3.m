function [u,s,v]=my_svd_3(A)
    [m,n]=size(A);
    [u,b,v,r]=bidiagonal_r_guess(A,n/2+1,m,n);
    [u(:,1:n),s,v]=dk_svd(b,u(:,1:n),v,n);
    [s,v]=change_signval(s,v,r);
    s=[s;zeros(m-n)];

end