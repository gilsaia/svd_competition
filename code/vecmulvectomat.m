function [S] = vecmulvectomat(x,y)
    [m,n]=size(x);
    [t,z]=size(y);
    S=zeros(m,z);
    for i=1:z
        b=y(i);
        for j=1:m
            S(j,i)=x(j)*b;
        end
    end 
end