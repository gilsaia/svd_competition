function [S,V] = change_signval(S,V,r)
    for i=1:r
        if S(i,i)<0
            S(i,i)=S(i,i)*-1;
            V(:,i)=scalemat(-1,V(:,i));
        end
    end
end