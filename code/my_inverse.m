function [inv_AA] = my_inverse(A)
    [m,n]=size(A);
    I=eye(m,m);
    [u,b,r]=bidiagonal_r_guess_withoutv(A,n/2+1,m,n);
    [u(:,1:n),s]=dk_svd_withoutv(b,u(:,1:n),n);
    ut=u';
    for i=1:r
        st=s(i,i)*s(i,i)';
        vt=-st/(1+st);
        u(:,i)=scalemat(vt,u(:,i));
    end
    u=matmulf(u(:,1:r),ut(1:r,:));
    inv_AA=u+I;
end