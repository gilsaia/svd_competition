function [w]=sv_approx(s,delta,r,n)
    w=zeros(n,1);
    presum=0;
    for i=1:n
        presum=presum+s(i);
    end
    delta(n+1)=delta(n)+presum;
    for i=n:-1:n-r+1
        y=get_init(i,s,delta,n);
        % disp(y);
        cnt=0;
        while 1
            cnt=cnt+1;
            eta=get_update(i,y,n,s,delta);
            if cnt>10000
                disp('error in iter')
                break
            end
            if ~stop(i,eta,delta,y)
                y=eta+y;
            else
                break;
            end
        end
        w(n-i+1)=y;
    end
end


function [is_stop]=stop(k,eta,delta,y)
    v=100*eps*min(abs(delta(k)-y),abs(delta(k+1)-y));
    % disp(abs(eta));
    % disp(v);
    if abs(eta)<=v
        is_stop=1;
    else
        is_stop=0;
    end
end


function [eta]=get_update(k,y,n,s,delta)
    if k<n
        [fy,~]=f(k,y,s,delta,n);
        [dfy,dpy,dqy]=df(k,y,s,delta,n);
        delta1=delta(k)-y;
        delta2=delta(k+1)-y;
        a=(delta1+delta2)*fy-delta1*delta2*dfy;
        b=delta1*delta2*fy;
        c=fy-delta1*dpy-delta2*dqy;
        if a<=0
            eta=(a-sqrt(a^2-4*b*c))/(2*c);
        else
            eta=2*b/(a+sqrt(a^2-4*b*c));
        end
    else
        [fy,~]=f(n-1,y,s,delta,n);
        [dfy,dpy,~]=df(n-1,y,s,delta,n);
        delta1=delta(n-1);
        delta2=delta(n);
        a=(delta1+delta2)*fy-delta1*delta2*dfy;
        b=delta1*delta2*fy;
        c=fy-delta1*dpy-s(n)/delta2;
        if a>=0
            eta=(a+sqrt(a^2-4*b*c))/(2*c);
        else
            eta=2*b/(a-sqrt(a^2-4*b*c));
        end
    end
end


function [y]=get_init(k,s,delta,n)
    tmp=(delta(k)+delta(k+1))/2;
    if k<n
        [F,G]=f(k,tmp,s,delta,n);
        if F>=0
            K=k;
            a=G*(delta(k+1)-delta(k))+s(k)+s(k+1);
            b=s(k)*(delta(k+1)-delta(k));
        else
            K=k+1;
            a=-G*(delta(k+1)-delta(k))+s(k)+s(k+1);
            b=-s(k+1)*(delta(k+1)-delta(k));
        end
        if a>0
            tau=2*b/(a+sqrt(a^2-4*b*G));
        else
            tau=(a-sqrt(a^2-4*b*G))/(2*G);
        end
    else
        K=n;
        [F,G]=f(n-1,tmp,s,delta,n);
        if F<=0 && 1+G<=-s(n-1)/(delta(n-1)-delta(n+1))-s(n)/(delta(n)-delta(n+1))
            tau=delta(n+1);
        else
            a=-G*(delta(n)-delta(n-1))+s(n-1)+s(n);
            b=-s(n)*(delta(n)-delta(n-1));
            if a>=0
                tau=(a+sqrt(a^2-4*b*G))/(2*G);
            else
                tau=2*b/(a-sqrt(a^2-4*b*G));
            end
        end
    end
    y=min(tau+delta(K),(delta(n)+delta(n+1))/2);
end


function [F,G]=f(k,y,s,delta,n)
    F=1;
    diff=0;
    for i=1:n
        F=F+s(i)/(delta(i)-y);
        if i==k || i==k+1
            diff=diff+s(i)/(delta(i)-y);
        end
    end
    G=F-diff;
end


function [dfy,dpy,dqy]=df(k,y,s,delta,n)
    dfy=0;
    dpy=0;
    for i=1:n
        dfy=dfy+s(i)/(delta(i)-y)^2;
        if i<=k
            dpy=dpy+s(i)/(delta(i)-y)^2;
        end
    end
    dqy=dfy-dpy;
end