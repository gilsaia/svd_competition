function [u,s,v]=my_svd_1(A,r)
    
    tic;
    [U,B,V]=bidiagonal_new(A,r);
    bitime=toc;
    [m,n]=size(B);
    Bo=B(1:n,:);
    disp(sprintf('Bidiagnoal time %f',bitime));
    tic;
    [Ux,Bt,Vx]=bidiagnoal_svd(Bo);
    bisvdtime=toc;
    disp(sprintf('Bidiagnoal time %f,SVD time %f', bitime,bisvdtime));
    remain=Bo-Ux'*Bt*Vx;
    e1=norm(remain,'fro')/norm(Bo,'fro');
    g1=max(norm(eye(n)-Ux*Ux','fro'),norm(eye(n)-Vx*Vx','fro'));
    disp(sprintf('E1:%f,G1:%f',e1,g1));
    [u,s,v]=svd(A);

end