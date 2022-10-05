function [C] = checkmat()
    A=rand(17,19);
    B=rand(19,23);
    tic;
    C=matmul(A,B);
    b=toc;
    c=norm(C-A*B,'fro');
    disp(sprintf('Time b:%f,Norm:%f', b,c));
end