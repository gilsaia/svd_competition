function [A] = generate_matrix()
    r1=rand(1,32);
    r2=zeros(1,96);
    r=[r1,r2];
    S1=diag(r);
    S2=zeros(128,128);
    S=[S1;S2];
    R=rand(256,128);
    [U,s,V]=svd(R);
    A=U*S*V;
end