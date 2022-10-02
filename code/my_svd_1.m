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


    e1=norm(U*B*V-A,'fro')/norm(A,'fro');
    disp(sprintf('bidiagonalization E1:%f',e1));

    remain=Bo-Ux'*Bt*Vx;
    e1=norm(remain,'fro')/norm(Bo,'fro');
    g1=max(norm(eye(n)-Ux*Ux','fro'),norm(eye(n)-Vx*Vx','fro'));
    disp(sprintf('bidiagonal svd E1:%f,G1:%f',e1,g1));


    A_rec = U(:,1:n)*Ux'*Bt*Vx*V;
    e1=norm(A_rec-A,'fro')/norm(A,'fro');
    disp(sprintf('all process E1:%f,G1:%f',e1,g1));


    [u,s,v] = svd(A);
    e1=norm(u*s*v'-A,'fro')/norm(A,'fro');
    disp(sprintf('build-in E1:%f,G1:%f',e1,g1));
end