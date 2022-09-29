function [u,s,v]=my_svd_1(A,r)
    
    [u,s,v]=jacobi_svd(A);
    % [u,s,v]=svd(A);

end