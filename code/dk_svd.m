function [U,S,V,B] = dk_svd(B,U,V,n)
    % B=U*S*V'
    % input param is U from bid with size(m,m) cut to size(m,n)
    theta=1e-19;
    tol=1e-18;
    upperd=n;
    while 1
        e=upperd;
        while e>1&&abs(B(e-1,e))<1e-15
            e-=1;
            upperd-=1;
        end
        s=e-1;
        while s>1&&abs(B(s-1,s))>1e-15
            s-=1;
        end
        if e==1
            break
        end
        Bt=B(s:e,s:e);
        [Bt,mumin]=stop_criterion(Bt,tol);
        if (e-s)==1
            [U,Bt,V]=twoeleQR(Bt,U,V,s);
        else
            bmax=max(max(Bt));
            upper=max(bmax/2,1e-16);
            lower=min(mumin*n^(1/2),mumin*n^-(1/2));
            if n*lower/upper<max(theta/tol,1e-1)
                [U,Bt,V]=implicitQR(Bt,U,V,s);
            else
                [U,Bt,V]=standardQR(Bt,U,V,s);
            end
        end
        % disp(sprintf('E:%d S:%d',e,s));
        B(s:e,s:e)=Bt;
    end
    S=diag(diag(B));
end

function [U,B,V] = twoeleQR(B,U,V,st)
    a=B(1,1);
    b=B(1,2);
    c=B(2,2);
    abcsum=a^2+b^2+c^2;
    delta=sqrt(abcsum^2-4*a^2*c^2);
    lambda=zeros(2,1);
    if(abs(abcsum+delta)>abs(abcsum-delta))
        lambda(1)=(abcsum+delta)/2;
        lambda(2)=(abcsum-delta)/2;
    else
        lambda(1)=(abcsum-delta)/2;
        lambda(2)=(abcsum+delta)/2;
    end

    v=ones(2,1);
    v(2)=-(a^2-lambda(1))/(a*b);

    vnorm=norm(v,'fro');
    v(1)/=vnorm;
    v(2)/=vnorm;

    u=ones(2,1);
    u(2)=-(a^2-lambda(2))/(a*b);

    unorm=norm(u,'fro');
    u(1)/=unorm;
    u(2)/=unorm;

    vsum=v(1)^2+v(2)^2;
    pot=u(1)*v(1)+u(2)*v(2);
    u(1)=u(1)-(pot/vsum)*v(1);
    u(2)=u(2)-(pot/vsum)*v(2);

    vt=[v,u];
    V(:,st:st+1)=matmul(V(:,st:st+1),vt);

    ut=matmul(B,vt);
    ut=matmul(ut,diag([1/sqrt(lambda(1)),1/sqrt(lambda(2))]));
    U(:,st:st+1)=matmul(U(:,st:st+1),ut);

    lambda(1)=sqrt(lambda(1));
    lambda(2)=sqrt(lambda(2));
    B=diag(lambda);
end

function [Bt,mumin] = stop_criterion(Bt,tol)
    [n,n]=size(Bt);
    mu=zeros(n,1);
    mu(1)=Bt(1,1);
    mumin=mu(1);
    for j=1:n-1
        if abs(Bt(j,j+1))<mu(j)*tol
            Bt(j,j+1)=0;
            return
        end
        mu(j+1)=abs(Bt(j+1,j+1))*mu(j)/(mu(j)+abs(Bt(j,j+1)));
        mumin=min(mumin,mu(j+1));
    end
end


function [cs,sn,r] = rot(f,g)
    if f==0
        cs=0;
        sn=1;
        r=g;
    elseif abs(f)>abs(g)
        t=g/f;
        tt=sqrt(1+t^2);
        cs=1/tt;
        sn=t*cs;
        r=f*tt;
    else
        t=f/g;
        tt=sqrt(1+t^2);
        sn=1/tt;
        cs=t*sn;
        r=g*tt;
    end
end

% B*[c -s;s c]
function [b] = updateBleft(b,c,s)
    n=size(b,1);
    for i=1:n
        tmp=b(i,1);
        b(i,1)=c*tmp+s*b(i,2);
        b(i,2)=-s*tmp+c*b(i,2);
    end
end

% [c s;-s c]*B
function [b] = updateBright(b,c,s)
    n=size(b,2);
    for i=1:n
        tmp=b(1,i);
        b(1,i)=c*tmp+s*b(2,i);
        b(2,i)=-s*tmp+c*b(2,i);
    end
end

function [U,B,V] = implicitQR(B,U,V,st)
    n=size(B,1);
    for i=1:n-1
        [c s r]=rot(B(i,i),B(i,i+1));
        B(max(i-1,1):i+1,i:i+1)=updateBleft(B(max(i-1,1):i+1,i:i+1),c,s);
        % B(:,i:i+1)=matmul(B(:,i:i+1),[c -s;s c]);
        V(:,st+i-1:st+i)=matmul(V(:,st+i-1:st+i),[c -s;s c]);
        [c s r]=rot(B(i,i),B(i+1,i));
        B(i:i+1,i:min(i+2,n))=updateBright(B(i:i+1,i:min(i+2,n)),c,s);
        % B(i:i+1,:)=matmul([c s;-s c],B(i:i+1,:));
        U(:,st+i-1:st+i)=matmul(U(:,st+i-1:st+i),[c -s;s c]);
    end
end


function [U,B,V] = standardQR(B,U,V,st)
    n=size(B,1);
    a=B(n-1,n-1);
    b=B(n,n);
    c=B(n-1,n);
    d=(a^2-b^2)/2;
    mu=b-c^2/(d+sign(d)*sqrt(d^2+c^2));
    x=B(1,1)^2-mu;
    z=B(1,1)*B(1,2);
    for i=1:n-1
        [c s r]=rot(x,z);
        B(max(i-1,1):i+1,i:i+1)=updateBleft(B(max(i-1,1):i+1,i:i+1),c,s);
        % B(:,i:i+1)=matmul(B(:,i:i+1),[c -s;s c]);
        V(:,st+i-1:st+i)=matmul(V(:,st+i-1:st+i),[c -s;s c]);
        x=B(i,i);
        z=B(i+1,i);
        [c s r]=rot(x,z);
        B(i:i+1,i:min(i+2,n))=updateBright(B(i:i+1,i:min(i+2,n)),c,s);
        % B(i:i+1,:)=matmul([c s;-s c],B(i:i+1,:));
        U(:,st+i-1:st+i)=matmul(U(:,st+i-1:st+i),[c -s;s c]);
        if i<n-1
            x=B(i,i+1);
            z=B(i,i+2);
        end
    end
end