function [u,s,v]=my_svd_1(A,r)

    [U,B,V]=checkbid(A);
    
    [u,s,v]=svd(A);

end