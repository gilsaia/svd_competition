function [e1,g1] = check_svd(A,U,S,V)
    AT=U*S*V';
    Remain=A-AT;
    e1=norm(Remain,'fro')/norm(A,'fro');
    ut=U'*U;
    eye_u=eye(size(ut));
    re_u=eye_u-ut;
    g1_u=norm(re_u,'fro')/norm(eye_u,'fro');
    vt=V'*V;
    eye_v=eye(size(vt));
    re_v=eye_v-vt;
    g1_v=norm(re_v,'fro')/norm(eye_v,'fro');
    g1=max(g1_u,g1_v);
    disp(sprintf('E1:%f,G1:%f',e1,g1));
end