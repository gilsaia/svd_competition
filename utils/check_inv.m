function [e2] = check_inv(A,inv)
    AT=A*A';
    [m,m]=size(AT);
    res=eye(m)+AT;
    e2=norm(res*inv-eye(m),'fro')/sqrt(m);
end