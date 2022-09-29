function [u,s,v]=my_svd_1(A,r)
    
    [u,sigma,v]=jacobi_svd2(A,'none');
    size(u)
    size(sigma)
    size(v)
    s=diag(sigma);
    % [m,n]=size(A);
    % temp=zeros(m-n,n);
    % s=[s;temp];
    % [u,s,v]=jacobi_svd(A);
    % [u,s,v]=svd(A);

end