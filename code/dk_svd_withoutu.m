function [S,V] = dk_svd_withoutu(B,V,n)
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
            [Bt,V]=twoelementQR(Bt,V,s);
        else
            bmax=max(max(Bt));
            upper=max(bmax/2,1e-16);
            lower=min(mumin*n^(1/2),mumin*n^-(1/2));
            if n*lower/upper<max(theta/tol,1e-5)
                [Bt,V]=implicitQR(Bt,V,s);
            else
                [Bt,V]=standardQR(Bt,V,s);
            end
        end
        B(s:e,s:e)=Bt;
    end
    S=diag(diag(B));
end

function [B,V] = twoelementQR(B,V,st)
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

function [B,V] = implicitQR(B,V,st)
    n=size(B,1);
    for i=1:n-1
        [c s r]=rot(B(i,i),B(i,i+1));
        B(max(i-1,1):i+1,i:i+1)=updateBleft(B(max(i-1,1):i+1,i:i+1),c,s);
        V(:,st+i-1:st+i)=matmulf(V(:,st+i-1:st+i),[c -s;s c]);
        [c s r]=rot(B(i,i),B(i+1,i));
        B(i:i+1,i:min(i+2,n))=updateBright(B(i:i+1,i:min(i+2,n)),c,s);
    end
end

function [B,V] = standardQR(B,V,st)
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
        if i<n-1
            x=B(i,i+1);
            z=B(i,i+2);
        end
    end
end