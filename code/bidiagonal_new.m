function [U,B,V] = bidiagonal_new(A,r)
    [m,n]=size(A);
    B=zeros(m,n);
    d=zeros(n,1);
    e=zeros(n,1);
    U=eye(m);
    V=eye(n);
    for j=1:n
        [alpha,tau,v]=householder_lapack(A(j:end,j),m-j+1);
        d(j)=real(alpha); 
        Q = eye(m);
        Q(j:end,j:end) = Q(j:end,j:end)-tau*v*v';
        U=U*Q;
        A=Q'*A;

        if j<n
            [alpha,tau,v]=householder_lapack(A(j,j+1:end)',n-j);
            e(j)=real(alpha);
            P=eye(n);
            P(j+1:end,j+1:end)=P(j+1:end,j+1:end)-tau*v*v';
            V=V*P;
            A=A*P;
        end
    end
    B=real(A);
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
    v(2:end)=A(2:end)*alpha;
    alpha=beta;
end

function [sigma,w] = householderp(x,len_x)
    x_norm = norm(x, 'fro');
    niu = sign(real(x(1))) * x_norm;
    k = (abs(real(x(1)))+x_norm) / x_norm;
    xt=x;
    xt(1)+=niu;
    w = xt;
    sigma = k / (niu*k*(x(1)+niu)');
end

function [delta] = householder(x, len_x)
    x_norm = norm(x, 'fro');
    niu = sign(real(x(1))) * x_norm;
    k = (abs(real(x(1)))+x_norm) / x_norm;
    e1 = zeros(len_x,1);
    e1(1) = 1;
    w = (x+niu*e1) * sqrt(k) / (x(1)+niu);
    sigma = (x(1)+niu) / (niu*k);
    delta = sigma * w * w';
end

function [sigma,w] = householderq(x,len_x)
    x_norm = norm(x, 'fro');
    niu = sign(real(x(1))) * x_norm;
    xt=x;
    xt(1)+=niu;
    w = xt;
    sigma = 1 / (niu*(x(1)+niu)');
end

function [delta] = householder2(x, len_x)
    x_norm = norm(x, 'fro');
    niu = sign(real(x(1))) * x_norm;
    e1 = zeros(len_x,1);
    e1(1) = 1;
    w = (x+niu*e1) / (x(1)+niu);
    sigma = (x(1)+niu) / niu;
    delta = sigma * w * w';
end
