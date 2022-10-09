function [U,S,V] = dc_svd(B,n,r)
    % B is down trian(B') and add zeros
    if r<2
        % need to recovery B
        Bt=B';
        Bt=[Bt;zeros(1,n+1)];
        [S,V]=dk_svd_withoutu(Bt,eye(n+1),n+1);
        [S,V]=change_signval(S,V,r);
        U=V;
        return
    end
    tk=round(r/2);
    [U1,S1,V1]=dc_svd(B(1:tk,1:tk-1),tk-1,tk-1);
    [U2,S2,V2]=dc_svd(B(tk+1:end,tk+1:end),n-tk,r-tk);
    alphak=B(tk,tk);
    betak=B(tk+1,tk);
    lambda1=U1(tk,tk);
    l1=U1(tk,1:tk-1)';
    lambda2=U2(1,n-tk+1);
    l2=U2(1,1:n-tk)';
    sigma=zeros(n,1);
    d=zeros(n,1);
    sigma(1)=sqrt((alphak*lambda1)^2+(betak*lambda2)^2);
    for i=1:tk-1
        sigma(i+1)=(alphak*l1(i));
    end
    for i=1:n-tk
        sigma(i+tk)=(betak*l2(i));
    end
    for i=1:tk-1
        d(i+1)=S1(i,i);
    end
    for i=1:n-tk
        d(i+tk)=S2(i,i);
    end
    [ds,dindex]=sort(d(2:n));
    d(2:n)=ds;
    sigmas=sigma;
    for i=2:n
        sigmas(i)=sigma(dindex(i-1)+1);
    end
    sigma2=sigmas;
    for i=1:n
        sigma2(i)=sigmas(i)^2;
    end
    d2=d;
    for i=1:n
        d2(i)=d(i)^2;
    end
    [w]=sv_approx(sigma2,d2,r,n);
    disp('Sigma');
    disp(sigmas);
    disp('D');
    disp(d);
    for i=r+1:n
        w(i)=0;
    end
    U=zeros(n,n);
    V=zeros(n,n);
    for i=1:n
        ui=zeros(n,1);
        for j=1:n
            ui(j)=sigmas(j)/(d2(j)-w(i));
        end
        uino=norm(ui,'fro');
        ui=scalemat(1/uino,ui);
        U(:,i)=ui;
    end
    for i=1:n
        vi=zeros(n,1);
        vi(1)=-1;
        for j=2:n
            vi(j)=d(j)*sigmas(j)/(d2(j)-w(i));
        end
        vino=norm(vi,'fro');
        vi=scalemat(1/vino,vi);
        V(:,i)=vi;
    end
    signval=zeros(n,1);
    for i=1:r
        signval(i)=sqrt(w(i));
    end
    U=[U zeros(n,1);zeros(1,n) 1];
    S=diag(signval);
end