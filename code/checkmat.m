function [C] = checkmat(A)
    tic;
    B=matmul(A,A');
    a=toc;
    tic;
    C=matmulf(A,A');
    b=toc;
    c=norm(B-C,'fro');
    disp(sprintf('Time a:%f,Time b:%f,Norm:%f', a,b,c));
end