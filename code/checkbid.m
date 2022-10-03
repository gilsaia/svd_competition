function [U,B,V] = checkbid(A,r)
    [U,B,V]=bidiagonal_new(A,r);
    AT=U*B*V';
    remain=A-AT;
    check=norm(remain,'fro')/norm(A,'fro');
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
    elesum=B(n,n);
    allsum=sum(B(n,1:n));
    disp(sprintf('index:%d,twoelesum:%f allsum:%f',n,elesum,allsum));
    if (allsum-elesum)>1e-5
        disp('Sum Wrong!');
        return
    end
    for i=n+1:m
        allsum=sum(B(i,1:n));
        disp(sprintf('index:%d,allsum:%f',i,allsum));
        if allsum>1e-5
            disp('Sum Wrong!');
            return
        end        
    disp('Bidiagnoal Pass!');
end