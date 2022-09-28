function [U,B,V] = bidiagnoal(A)
    [m,n]=size(A);
    B=A;
    U=eye(m);
    V=eye(n);
    for i=1:n-1
        x=householder(B(:,i),i,m);

        % bt=x'*B;
        % xt=2*x;
        % B=B-xt*bt;
        bt=matmul(x',B);
        xt=scalemat(2,x);
        B=B-vecmulvectomat(xt,bt);

        % ut=U*xt;
        % U=U-ut*x';
        ut=matmul(U,xt);
        U=U-vecmulvectomat(ut,x');

        y=householder(B(i,:)',i+1,n);

        % yt=2*y;
        % bt=B*yt;
        % B=B-bt*y';
        yt=scalemat(2,y);
        bt=matmul(B,yt);
        B=B-vecmulvectomat(bt,y');

        % vt=y'*V;
        % V=V-yt*vt;
        vt=matmul(y',V);
        V=V-vecmulvectomat(yt,vt);
    end
    x=householder(B(:,n),n,m);

    % bt=x'*B;
    % xt=2*x;
    % B=B-xt*bt;
    bt=matmul(x',B);
    xt=scalemat(2,x);
    B=B-matmul(xt,bt);

    % ut=U*xt;
    % U=U-ut*x';
    ut=matmul(U,xt);
    U=U-matmul(ut,x');
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

% function [Q] = findQ(B,k,m)
%     x=zeros(m,1);
%     bkk=B(k,k);
%     bkka=abs(bkk);
%     sk=norm(B(k:m,k),'fro');
%     x(k)=sqrt((1+bkka/sk)/2);
%     ck=(2*sk*bkk/bkka*x(k))^-1;
%     for i=k+1:m
%         x(i)=ck*B(i,k);
%     end
%     Q=eye(m);
%     Q=Q-2*(x*x');
%     % Q=Q-scalemat(2,vecmulvectomat(x));
% end

% function [P] = findP(B,k,n)
%     y=zeros(n,1);
%     bkk1=B(k,k+1)';
%     bkk1a=abs(bkk1);
%     tk=norm(B(k,k+1:n)','fro');
%     y(k+1)=sqrt((1+bkk1a/tk)/2);
%     dk=(2*tk*bkk1/bkk1a*y(k+1))^-1;
%     for i=k+2:n
%         y(i)=dk*conj(B(k,i)');
%     end
%     P=eye(n);
%     P=P-2*(y*y');
%     % P=P-scalemat(2,vecmulvectomat(y));
% end