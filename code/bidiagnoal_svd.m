function [U,S,V] = bidiagnoal_svd(B)
    % B=U*S*V'
    [n,n]=size(B);
    theta=1e-16;
    tol=1e-14;
    U=eye(n);
    V=eye(n);
    while 1
        s=1;
        e=n;
        while s<n&&abs(B(s,s+1))<1e-8
            s+=1;
        end
        while e>1&&abs(B(e-1,e))<1e-8
            e-=1;
        end
        if s>=e
            break
        end
        Bt=B(s:e,s:e);
        [Bt,mumin]=stop_criterion(Bt,tol);
        bmax=max(max(Bt));
        upper=bmax/2;
        lower=min(mumin*n^(1/2),mumin*n^-(1/2));
        if n*lower/upper<max(theta/tol,1e-2)
            [U,Bt,V]=implicitQR(Bt,U,V,s);
        else
            [U,Bt,V]=standardQR(Bt,U,V,s);
        end
        B(s:e,s:e)=Bt;
    end
    S=diag(diag(B));
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


function [U,B,V] = implicitQR(B,U,V,st)
    n=size(B,1);
    for i=1:n-1
        [c s r]=rot(B(i,i),B(i,i+1));
        Q=eye(n);
        Q(i:i+1,i:i+1)=[c s;-s c];
        B(:,i:i+1)=matmul(B(:,i:i+1),[c -s;s c]);
        V(:,st+i-1:st+i)=matmul(V(:,st+i-1:st+i),[c -s;s c]);
        [c s r]=rot(B(i,i),B(i+1,i));
        Q=eye(n);
        Q(i:i+1,i:i+1)=[c s;-s c];
        B(i:i+1,:)=matmul([c s;-s c],B(i:i+1,:));
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
        Q=eye(n);
        Q(i:i+1,i:i+1)=[c s;-s c];
        B(:,i:i+1)=matmul(B(:,i:i+1),[c -s;s c]);
        V(:,st+i-1:st+i)=matmul(V(:,st+i-1:st+i),[c -s;s c]);
        x=B(i,i);
        z=B(i+1,i);
        [c s r]=rot(x,z);
        Q=eye(n);
        Q(i:i+1,i:i+1)=[c s;-s c];
        B(i:i+1,:)=matmul([c s;-s c],B(i:i+1,:));
        U(:,st+i-1:st+i)=matmul(U(:,st+i-1:st+i),[c -s;s c]);
        if i<n-1
            x=B(i,i+1);
            z=B(i,i+2);
        end
    end
end