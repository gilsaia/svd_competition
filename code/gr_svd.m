function [U,S,V] = gr_svd(B)
    [n,n]=size(B);
    sigma=1e-14
    while 1
        for i=1:n-1
            if abs(B(i,i+1))<=sigma*(abs(B(i,i))+abs(B(i+1,i+1)))
                B(i,i+1)=0;
            end
        end
        e=n;
        s=0;
        while e>0 && B(e-1,e)==0
            e-=1;
        end
        while s<n && B(s,s+1)==0
            s+=1;
        end
        if s>=e
            break
        end
        flag=0;
        e=e-1;
        for i=s:e
            if B(i,i)==0
                B(i,i+1)=0;
                flag=1;
            end
        end
        if flag==0
            [B,U,V]=gr_step(B,U,V,s,e);
        end
    end
    S=B;
end

function [B,Q,P] = gr_step(B,Q,P,s,e)
    a=B(e-1,e-1);
    b=B(e-1,e);
    c=B(e,e);
    lamt=(a^2+b^2+c^2);
    delta=lamt^2-4*a^2*c^2;
    lam1=(lamt+delta)/2;
    lam2=(lamt-delta)/2;
    comp=b^2+c^2;
    remain1=abs(lam1-comp);
    remain2=abs(lam2-comp);
    sigma=B(s,s)*B(s,s+1);
    if remain1>remain2
        y=B(s,s)^2-remain2;
    else
        y=B(s,s)^2-remain1;
    end
    for k=s:e
        tmpsum=sqrt(y^2+sigma^2);
        c=y/tmpsum;
        s=-sigma/tmpsum;
        B=B;
        P=P;
        y=B(k,k);
        sigma=B(k+1,k);
        tmpsum=sqrt(y^2+sigma^2);
        c=y/tmpsum;
        s=-sigma/tmpsum;
        B=B;
        Q=Q;
        if k<e-1
            y=B(k,k+1);
            sigma=B(k,k+2);
        end
    end
end