function [C] = matmulf(A,B)
    [m,t]=size(A);
    [t,n]=size(B);
    C=zeros(m,n);
    for i=1:m
        for k=1:t
            s=A(i,k);
            for j=1:n
                C(i,j)+=s*B(k,j);
            end
        end
    end
end