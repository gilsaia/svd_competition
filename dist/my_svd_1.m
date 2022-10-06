function [u,s,v]=my_svd_1(A,r)
    [m,n]=size(A);
    if m <= 256
        bound=50;
    elseif m<=512
        bound=40;
    else
        bound=30;
    end
    [u,b,v]=bidiagonal_r(A,min(r,bound),m,n);
    [u(:,1:n),s,v]=dk_svd(b,u(:,1:n),v,n);
    [s,v]=change_signval(s,v,r);
    s=[s;zeros(m-n)];
end