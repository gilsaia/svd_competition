function [U,S,V] = dc_svd(B,n,r)
    % B is down trian(B') and add zeros
    if r<256
        % need to recovery B
        Bt=B';
        Bt=Bt(1:n,1:n);
        [U,S,V]=dk_svd(Bt,eye(n),eye(n),n);
        return
    end
    tk=round(r/2);
    [U1,S1,V1]=dc_svd(B(1:tk,1:tk-1),tk-1,tk-1);
    [U2,S2,V2]=dc_svd(B(tk+1:end,tk+1:end),n-tk,r-tk-1);
    alphak=B(tk,tk);
    betak=B(tk+1,tk);
    lambda1=U1(tk,tk);
    l1=U1(tk,1:tk-1)';
    lambda2=U2(1,n-tk+1);
    l2=U2(1,1:n-tk)';
    sigma=zeros(n,1);
    d=zeros(n,1);
    sigma(1)=sqrt((alphak*lambda1)^2+(betak*lambda2)^2);
    d(1)=sigma(1);
    for i=1:tk-1
        sigma(i+1)=alphak*l1(i);
    end
    for i=1:n-tk
        sigma(i+tk)=beta*l2(i);
    end
    for i=1:tk-1
        d(i+1)=S1(i,i);
    end
    for i=1:n-tk
        d(i+tk)=S2(i,i);
    end
    
    
end