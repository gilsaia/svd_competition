function [S] = vecmulvectomat(x,y)
    [m,n]=size(x);
    [t,z]=size(y);
    S=zeros(m,z);
    for i=1:m
        for j=1:z
            S(i,j)=x(i)*y(j);
        end
    end 
end