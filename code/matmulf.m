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
                y+=1;
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
                m+=1;
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
            m+=1;
        end
        if mod(n,2)==1
            A=[A zeros(m,1)];
            n+=1;
            B=[B;zeros(1,y)];
            x+=1;
        end
        if mod(y,2)==1
            B=[B zeros(x,1)];
            y+=1;
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
                c+=AT(k,i)*B(k,j);
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