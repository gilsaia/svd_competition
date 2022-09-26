function [U,B,V] = checkbid(A)
    [U,B,V]=bidiagnoal(A);
    AT=U*B*V';
    remain=A-AT;
    check=norm(remain,'fro');
    disp(sprintf('norm:%f',check));
    if check>1e-5
        disp('Norm Wrong!');
        % return
    end
    [m,n]=size(B);
    for i=1:n-1
        elesum=B(i,i)+B(i,i+1);
        allsum=sum(B(i,1:n));
        disp(sprintf('index:%d,twoelesum:%f allsum:%f',i,elesum,allsum));
        if (allsum-elesum)>1e-5
            disp('Sum Wrong!');
            return
        end
    end
    disp('Bidiagnoal Pass!');
end