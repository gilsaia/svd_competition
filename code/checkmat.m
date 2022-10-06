function [C] = checkmat()
    A=rand(512,512);
    B=rand(512,512);
    tic;
    C=matmulf(A,B);
    b=toc;
    c=norm(C-A*B,'fro');
    disp(sprintf('Time b:%f,Norm:%f', b,c));
end