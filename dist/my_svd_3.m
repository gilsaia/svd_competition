function [u,s,v]=my_svd_3(A)
    [m,n]=size(A);
    [u,b,v,r]=bidiagonal_r_guess(A,n/2+1,m,n);
    [u(:,1:n),s,v]=dk_svd(b,u(:,1:n),v,n);
    [s,v]=change_signval(s,v,r);
    s=[s;zeros(m-n)];
end


function [U,B,V,r] = bidiagonal_r_guess(A,bound,m,n)
    B=zeros(m,n);
    d=zeros(n,1);
    e=zeros(n-1,1);
    U=eye(m);
    V=eye(n);
    r=bound;
    for j=1:bound
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
        if (abs(d(j))+abs(e(j)))<1e-10
            r=j;
            break
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
    tau=tau*(alpha*alpha');
    alpha=beta;
end


function [U,S,V] = dk_svd(B,U,V,n)
    % B=U*S*V'
    % input param is U from bid with size(m,m) cut to size(m,n)
    theta=1e-16;
    tol=1e-10;
    upperd=n;
    while 1
        e=upperd;
        while e>1&&abs(B(e-1,e))<1e-15
            e=e-1;
            upperd=upperd-1;
        end
        s=e-1;
        while s>1&&abs(B(s-1,s))>1e-15
            s=s-1;
        end
        if e==1
            break
        end
        Bt=B(s:e,s:e);
        [Bt,mumin]=stop_criterion(Bt,tol);
        if (e-s)==1
            [U,Bt,V]=twoelementQR(Bt,U,V,s);
        else
            bmax=max(max(Bt));
            upper=max(bmax/2,1e-16);
            lower=min(mumin*n^(1/2),mumin*n^-(1/2));
            if n*lower/upper<max(theta/tol,1e-5)
                [U,Bt,V]=implicitQR(Bt,U,V,s);
            else
                [U,Bt,V]=standardQR(Bt,U,V,s);
            end
        end
        B(s:e,s:e)=Bt;
    end
    S=diag(diag(B));
end

function [U,B,V] = twoelementQR(B,U,V,st)
    f=B(1,1);
    h=B(2,2);
    g=B(1,2);
    ft=B(1,1);
    fa=abs(ft);
    ht=B(2,2);
    ha=abs(ht);
    pmax=1;
    swap=0;
    if ha>fa
        swap=1;
        pmax=3;
        temp=ft;
        ft=ht;
        ht=temp;
        temp=fa;
        fa=ha;
        ha=temp;
    end
    gt=B(1,2);
    ga=abs(gt);
    if ga<1e-18
        ssmin=ha;
        ssmax=fa;
        clt=1;
        crt=1;
        slt=0;
        srt=0;
    else
        gasmal=1;
        if ga>fa
            pmax=2;
            if (fa/ga)<1e-16
                gasmal=0;
                ssmax=ga;
                if ha>1
                    ssmin=fa/(ga/ha);
                else
                    ssmin=(fa/ga)*ha;
                end
                clt=1;
                slt=ht/gt;
                srt=1;
                crt=ft/gt;
            end
        end
        if gasmal~=0
            d=fa-ha;
            if abs(d-fa)<1e-18
                l=1;
            else
                l=d/fa;
            end
            m=gt/ft;
            t=2-l;
            mm=m*m;
            tt=t*t;
            s=sqrt(tt+mm);
            if abs(l)<1e-18
                r=abs(m);
            else
                r=sqrt(l*l+mm);
            end
            a=0.5*(s+r);
            ssmin=ha/a;
            ssmax=fa*a;
            if abs(mm)<1e-18
                if abs(l)<1e-18
                    t=2*sign(ft)*sign(gt);
                else
                    t=gt/(d*sign(ft))+m/t;
                end
            else
                t=(m/(s+t)+m/(r+l))*(1+a);
            end
            l=sqrt(t*t+4);
            crt=2/l;
            srt=t/l;
            clt=(crt+srt*m)/a;
            slt=(ht/ft)*srt/a;
        end
    end
    if swap==1
        csl=srt;
        snl=crt;
        csr=slt;
        snr=clt;
    else
        csl=clt;
        snl=slt;
        csr=crt;
        snr=srt;
    end
    if pmax==1
        tsign=sign(csr)*sign(csl)*sign(f);
    elseif pmax==2
        tsign=sign(snr)*sign(csl)*sign(g);
    else
        tsign=sign(snr)*sign(snl)*sign(h);
    end
    ssmax=ssmax*sign(tsign);
    ssmin=ssmin*sign(tsign*sign(f)*sign(h));
    B=diag([ssmax,ssmin]);
    vt=[csr -snr;snr csr];
    V(:,st:st+1)=matmulf(V(:,st:st+1),vt);
    ut=[csl -snl;snl csl];
    U(:,st:st+1)=matmulf(U(:,st:st+1),ut);
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
        V(:,st+i-1:st+i)=matmulf(V(:,st+i-1:st+i),[c -s;s c]);
        [c s r]=rot(B(i,i),B(i+1,i));
        B(i:i+1,i:min(i+2,n))=updateBright(B(i:i+1,i:min(i+2,n)),c,s);
        U(:,st+i-1:st+i)=matmulf(U(:,st+i-1:st+i),[c -s;s c]);
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
        V(:,st+i-1:st+i)=matmulf(V(:,st+i-1:st+i),[c -s;s c]);
        x=B(i,i);
        z=B(i+1,i);
        [c s r]=rot(x,z);
        B(i:i+1,i:min(i+2,n))=updateBright(B(i:i+1,i:min(i+2,n)),c,s);
        U(:,st+i-1:st+i)=matmulf(U(:,st+i-1:st+i),[c -s;s c]);
        if i<n-1
            x=B(i,i+1);
            z=B(i,i+2);
        end
    end
end


function [S,V] = change_signval(S,V,r)
    for i=1:r
        if S(i,i)<0
            S(i,i)=S(i,i)*-1;
            V(:,i)=scalemat(-1,V(:,i));
        end
    end
end


function [X] = scalemat(t,X)
    [m,n]=size(X);
    for i=1:m
        for j=1:n
            X(i,j)=X(i,j)*t;
        end
    end
end


function [S] = vecmulvectomat(x,y)
    [m,n]=size(x);
    [t,z]=size(y);
    S=zeros(m,z);
    for i=1:z
        b=y(i);
        for j=1:m
            S(j,i)=x(j)*b;
        end
    end 
end


function [C] = matmulf(A,B)
    [m,n]=size(A);
    [x,y]=size(B);
    if n~=x
        error('err dim',A,B);
    end
    if m==2&&n==2
        if x>16||y>16
            orix=x;
            oriy=y;
            if mod(y,2)==1
                B=[B zeros(x,1)];
                y=y+1;
            end
            C=twotwoleft(A,B,x,y);
            C=C(1:orix,1:oriy);
            return
        end
    end
    if x==2&&y==2
        if m>16||n>16
            orim=m;
            orin=n;
            if mod(m,2)==1
                A=[A;zeros(1,n)];
                m=m+1;
            end
            C=twotwoleft(B',A',n,m);
            C=C(1:orin,1:orim)';
            return
        end
    end
    if m>=16&&n>=16&&y>=16
        orim=m;
        orin=n;
        orix=x;
        oriy=y;
        if mod(m,2)==1
            A=[A;zeros(1,n)];
            m=m+1;
        end
        if mod(n,2)==1
            A=[A zeros(m,1)];
            n=n+1;
            B=[B;zeros(1,y)];
            x=x+1;
        end
        if mod(y,2)==1
            B=[B zeros(x,1)];
            y=y+1;
        end
        C=strassen(A,B,m,n,y);
        C=C(1:orim,1:oriy);
        return
    end
    C=zeros(m,y);
    AT=A.';
    for i=1:m
        for j=1:y
            c=0;
            for k=1:n
                c=c+AT(k,i)*B(k,j);
            end
            C(i,j)=c;
        end
    end
end

function [C] = strassen(A,B,m,n,y)
    splitm=m/2;
    splitn=n/2;
    splity=y/2;
    A11=A(1:splitm,1:splitn);
    A12=A(1:splitm,splitn+1:end);
    A21=A(splitm+1:end,1:splitn);
    A22=A(splitm+1:end,splitn+1:end);
    B11=B(1:splitn,1:splity);
    B12=B(1:splitn,splity+1:end);
    B21=B(splitn+1:end,1:splity);
    B22=B(splitn+1:end,splity+1:end);
    M1=matmulf(A11+A22,B11+B22);
    M2=matmulf(A21+A22,B11);
    M3=matmulf(A11,B12-B22);
    M4=matmulf(A22,B21-B11);
    M5=matmulf(A11+A12,B22);
    M6=matmulf(A21-A11,B11+B12);
    M7=matmulf(A12-A22,B21+B22);
    C=[M1+M4-M5+M7 M3+M5;M2+M4 M1-M2+M3+M6];
end

function [C] = twotwoleft(A,B,x,y)
    a=A(1,1);
    b=A(1,2);
    c=A(2,1);
    d=A(2,2);
    splitx=x/2;
    splity=y/2;
    ta=scalemat(a,B(1:splitx,1:splity));
    tb=scalemat(b,B(splitx+1:end,1:splity));
    td=scalemat(a+b-c-d,B(splitx+1:end,splity+1:end));
    cd=B(1:splitx,splity+1:end)-B(splitx+1:end,splity+1:end);
    u=scalemat(c-a,cd);
    ca=B(1:splitx,splity+1:end)-B(1:splitx,1:splity);
    v=scalemat(c+d,ca);
    adc=B(1:splitx,1:splity)+B(splitx+1:end,splity+1:end)-B(1:splitx,splity+1:end);
    tadc=scalemat(c+d-a,adc);
    w=ta+tadc;
    bcad=B(splitx+1:end,1:splity)+B(1:splitx,splity+1:end)-B(1:splitx,1:splity)-B(splitx+1:end,splity+1:end);
    tbcad=scalemat(d,bcad);
    C=[ta+tb w+v+td;w+u+tbcad w+u+v];
end