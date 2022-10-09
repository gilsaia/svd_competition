function [u,s,v]=my_svd_3(A)
    [m,n]=size(A);
    [u,b,v,r]=bidiagonal_r_guess_energy(A,n/2+1,m,n);
    if r>256
        % prepare b
        bt=[b';zeros(1,n)];
        [ut,s,vt]=dc_svd(bt,n,r);
        % dc output ut is v,vt is u
        ut=ut(1:n,1:n);
        v=matmulf(v,ut);
        u=matmulf(u,[vt zeros(n,m-n);zeros(m-n,n) eye(m-n,m-n)]);
    else
        [u(:,1:n),s,v]=dk_svd(b,u(:,1:n),v,n);
    end
    [s,v]=change_signval(s,v,r);
    s=[s;zeros(m-n)];
end