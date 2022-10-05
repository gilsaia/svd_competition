function [U,B,V,d,e] = bidiagonal_r(A,r,m,n)
    B=zeros(m,n);
    d=zeros(n,1);
    e=zeros(n-1,1);
    U=eye(m);
    V=eye(n);
    for j=1:r
        [alpha,tau,v]=householder_lapack(A(j:end,j),m-j+1);
        d(j)=real(alpha); 

        tt=scalemat(tau,v);
        ut=matmulf(U(:,j:end),tt);
        U(:,j:end)=U(:,j:end)-vecmulvectomat(ut,v');

        bt=matmulf(A(j:end,j+1:end)',tt);
        A(j:end,j+1:end)=(A(j:end,j+1:end)'-vecmulvectomat(bt,v'))';

        if j<n
            [alpha,tau,v]=householder_lapack(A(j,j+1:end)',n-j);
            e(j)=real(alpha);

            tt=scalemat(tau,v);
            vt=matmulf(V(:,j+1:end),tt);
            V(:,j+1:end)=V(:,j+1:end)-vecmulvectomat(vt,v');

            bt=matmulf(A(j+1:end,j+1:end),tt);
            A(j+1:end,j+1:end)=A(j+1:end,j+1:end)-vecmulvectomat(bt,v');
        end
    end
    B=diag(d)+diag(e,1);
end

function [alpha,tau,v] = householder_lapack(A,n)
    alpha=A(1);
    xnorm=norm(A(2:end),'fro');
    alphr=real(alpha);
    alphi=imag(alpha);
    v=zeros(n,1);
    v(1)=1;
    if xnorm==0 && alphi==0
        tau=0;
        return
    end
    beta=-1*sign(alphr)*sqrt(alphr^2+alphi^2+xnorm^2);
    tau=complex((beta-alphr)/beta,-alphi/beta);
    alpha=1/(alpha-beta);
    v(1)=1/alpha;
    v(2:end)=A(2:end);
    tau*=(alpha*alpha');
    alpha=beta;
end