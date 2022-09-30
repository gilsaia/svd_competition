function [u,s,v]=my_svd_1(A,r)
    
    tic;
    [U,B,V]=bidiagnoal_r(A,r);
    bitime=toc;
    [m,n]=size(B);
    Bo=B(1:n,:);
    disp(sprintf('Bidiagnoal time %f',bitime));
    tic;
    [Ux,Bt,Vx]=bidiagnoal_svd(Bo);
    bisvdtime=toc;
    disp(sprintf('Bidiagnoal time %f,SVD time %f', bitime,bisvdtime));
    [u,s,v]=svd(A);

end