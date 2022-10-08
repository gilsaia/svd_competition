function [S,V] = change_signval(S,V,r)
    for i=1:r
        if S(r,r)<0
            S(r,r)=S(r,r)*-1;
            V(:,r)=scalemat(-1,V(:,r));
        end
    end
end