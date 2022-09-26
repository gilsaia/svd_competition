function [A] = scalemat(t,X)
    [m,n]=size(X);
    for i=1:m
        for j=1:n
            X(i,j)*=t;
        end
    end
    A=X;
end