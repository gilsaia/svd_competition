function [C] = checkmat()
    A=rand(1024,1024);
    B=rand(1024,1024);
    tic;
    C=matmulf(A,B);
    b=toc;
    c=norm(C-A*B,'fro');
    disp(sprintf('Time b:%f,Norm:%f', b,c));
end