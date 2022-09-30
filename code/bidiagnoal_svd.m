function [U,S,V] = bidiagnoal_svd(B)
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
            [U,Bt,V]=implicitQR(Bt,U,V);
        else
            [U,Bt,V]=implicitQR(Bt,U,V);
        end
        B(s:e,s:e)=Bt;
    end
    S=zeros(n);
    for i=1:n
        S(i,i)=B(i,i);
    end
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

function [v1c,v2c] = update(cs,sn,v1,v2)
    v1c=v1;
    v2c=v2;
    [t,n]=size(v1);
    for i=1:n
        v1c(i)=cs*v1(i)+sn*v2(i);
        v2c(i)=-sn*v1(i)+cs*v2(i);
    end
end

function [cs,sn,r] = ROT(f,g)
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

function [U,B,V] = implicitQR(B,U,V)
    [n,n]=size(B);
    oldcs=1;
    % oldsn=0;
    f=B(1,1);
    g=B(1,2);
    for i=1:n-1
        [cs,sn,r]=ROT(f,g);
        [V(i,:),V(i+1,:)]=update(cs,sn,V(i,:),V(i+1,:));
        if i~=1
            B(i-1,i)=oldsn*r;
        end
        f=oldcs*r;
        g=B(i+1,i+1)*sn;
        h=B(i+1,i+1)*cs;
        [cs,sn,r]=ROT(f,g);
        [U(i,:),U(i+1,:)]=update(cs,sn,U(i,:),U(i+1,:));
        B(i,i)=r;
        f=h;
        if i<n-1
            g=B(i+1,i+2);
        end
        oldcs=cs;
        oldsn=sn;
    end
    B(n-1,n)=h*sn;
    B(n,n)=h*cs;
end