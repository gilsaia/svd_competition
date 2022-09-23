function [e1,g1] = check_svd(A,U,S,V)
    AT=U*S*V';
    Remain=A-AT;
    e1=norm(Remain,'fro')/norm(A,'fro');
    eye_u=eye(size(U));
    re_u=eye_u-U'*U;
    g1_u=norm(re_u,'fro')/norm(eye_u,'fro');
    eye_v=eye(size(V));
    re_v=eye_v-V'*V;
    g1_v=norm(re_v,'fro')/norm(eye_v,'fro');
    g1=max(g1_u,g1_v);
end