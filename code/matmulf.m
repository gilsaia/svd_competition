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
    % BT=B';
    % tm=floor(m/4);
    % tn=floor(n/4);
    % ty=floor(y/4);
    % for ti=1:tm
    %     for tj=1:ty
    %         c=zeros(4,4);
    %         for tk=1:tn
    %             for i=1:4
    %                 al=A(ti*4-3:ti*4,tk*4-4+i);
    %                 a1=al(1);
    %                 a2=al(2);
    %                 a3=al(3);
    %                 a4=al(4);

    %                 bl=BT(tj*4-3:tj*4,tk*4-4+i);
    %                 c(1,1)+=a1*bl(1);
    %                 c(1,2)+=a1*bl(2);
    %                 c(1,3)+=a1*bl(3);
    %                 c(1,4)+=a1*bl(4);

    %                 c(2,1)+=a2*bl(1);
    %                 c(2,2)+=a2*bl(2);
    %                 c(2,3)+=a2*bl(3);
    %                 c(2,4)+=a2*bl(4);
                    
    %                 c(3,1)+=a3*bl(1);
    %                 c(3,2)+=a3*bl(2);
    %                 c(3,3)+=a3*bl(3);
    %                 c(3,4)+=a3*bl(4);

    %                 c(4,1)+=a4*bl(1);
    %                 c(4,2)+=a4*bl(2);
    %                 c(4,3)+=a4*bl(3);
    %                 c(4,4)+=a4*bl(4);
    %             end
    %         end
    %         for k=tn*4+1:n
    %             al=A(ti*4-3:ti*4,k);
    %             a1=al(1);
    %             a2=al(2);
    %             a3=al(3);
    %             a4=al(4);

    %             bl=BT(tj*4-3:tj*4,k);
    %             c(1,1)+=a1*bl(1);
    %             c(1,2)+=a1*bl(2);
    %             c(1,3)+=a1*bl(3);
    %             c(1,4)+=a1*bl(4);

    %             c(2,1)+=a2*bl(1);
    %             c(2,2)+=a2*bl(2);
    %             c(2,3)+=a2*bl(3);
    %             c(2,4)+=a2*bl(4);
                
    %             c(3,1)+=a3*bl(1);
    %             c(3,2)+=a3*bl(2);
    %             c(3,3)+=a3*bl(3);
    %             c(3,4)+=a3*bl(4);

    %             c(4,1)+=a4*bl(1);
    %             c(4,2)+=a4*bl(2);
    %             c(4,3)+=a4*bl(3);
    %             c(4,4)+=a4*bl(4);
    %         end
    %         C(ti*4-3:ti*4,tj*4-3:tj*4)=c;
    %     end
    % end
    % if mod(m,4)==0&&mod(y,4)==0
    %     return
    % end
    % for i=tm*4+1:m
    %     cl=zeros(1,y-ty*4);
    %     for tk=1:tn
    %         al=A(i,tk*4-3:tk);
    %         a1=al(1);
    %         a2=al(2);
    %         a3=al(3);
    %         a4=al(4);
    %         for j=ty*4+1:y
    %             bl=B(tk*4-3:tk,j);
    %             c=a1*bl(1)+a2*bl(2)+a3*bl(3)+a4*bl(4);
    %             cl(j-ty*4)+=c;
    %         end
    %     end
    %     for k=tn*4+1:n
    %         a=A(i,k);
    %         for j=ty*4+1:y
    %             cl(j-ty*4)+=a*BT(j,k);
    %         end
    %     end
    %     C(i,y-ty*4:y)=cl;
    % end

                
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