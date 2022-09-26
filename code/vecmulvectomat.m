function [S] = vecmulvectomat(x)
    [m,n]=size(x);
    S=zeros(m);
    for i=1:m
        for j=1:m
            S(i,j)=x(i)*x(j);
        end
    end 
end