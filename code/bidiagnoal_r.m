function [U,B,V] = bidiagnoal_r(A,r)
    [m,n]=size(A);
    B=A;
    U=eye(m);
    V=eye(n);
    for i=1:r
        x=householder(B(:,i),i,m);

        bt=x'*B;
        xt=2*x;
        B=B-xt*bt;
        % bt=matmul(x',B);
        % xt=scalemat(2,x);
        % B=B-vecmulvectomat(xt,bt);

        ut=U*xt;
        U=U-ut*x';
        % ut=matmul(U,xt);
        % U=U-vecmulvectomat(ut,x');

        y=householder(B(i,:)',i+1,n);

        yt=2*y;
        bt=B*yt;
        B=B-bt*y';
        % yt=scalemat(2,y);
        % bt=matmul(B,yt);
        % B=B-vecmulvectomat(bt,y');

        vt=y'*V;
        V=V-yt*vt;
        % vt=matmul(y',V);
        % V=V-vecmulvectomat(yt,vt);
    end
end

function [x] = householder(b,k,m)
    x=zeros(m,1);
    bkk=b(k);
    bkka=abs(bkk);
    sk=norm(b(k:m),'fro');
    x(k)=sqrt((1+bkka/sk)/2);
    ck=(2*sk*bkk/bkka*x(k))^-1;
    for i=k+1:m
        x(i)=ck*b(i);
    end
end