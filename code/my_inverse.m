function [inv_AA] = my_inverse(A)
    [m,n]=size(A);
    I=eye(m,m);
    [u,b,r]=bidiagonal_r_guess_withoutv(A,n/2+1,m,n);
    if r>256
        % prepare b
        bt=[b';zeros(1,n)];
        [s,vt]=dc_svd_withoutu(bt,n,r);
        % dc output ut is v,vt is u
        u=matmulf(u,[vt zeros(n,m-n);zeros(m-n,n) eye(m-n,m-n)]);
    else
        [u(:,1:n),s]=dk_svd_withoutv(b,u(:,1:n),n);
    end
    ut=u';
    for i=1:r
        st=s(i,i)*s(i,i)';
        vt=-st/(1+st);
        u(:,i)=scalemat(vt,u(:,i));
    end
    u=matmulf(u(:,1:r),ut(1:r,:));
    inv_AA=u+I;
end