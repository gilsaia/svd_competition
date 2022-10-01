function [U,B,V] = bidiagnoal_new(A,r)
    [m,n]=size(A);
    B=A;
    U=eye(m);
    V=eye(n);
    for j=1:n-1
        disp(sprintf('iter:%f',j));
        delta = householder2(B(j:end,j),m-j+1);        
        B(j:end,j:end) = B(j:end,j:end) - delta' * B(j:end,j:end);
        Q = eye(m);
        Q(j:end,j:end) = Q(j:end,j:end) - delta';
        U = matmul(U,Q);

        if j <= n-2
            delta = householder(B(j,j+1:end)',n-j);
            B(j:end,j+1:end) = B(j:end,j+1:end) - B(j:end,j+1:end) * delta;
            P = eye(n);
            P(j+1:end,j+1:end) = P(j+1:end,j+1:end) - delta;
            V = matmul(P,V);
        end
    end

end

function [delta] = householder(x, len_x)
    x_norm = norm(x, 'fro');
    niu = sign(real(x(1))) * x_norm;
    k = (abs(real(x(1)))+x_norm) / x_norm;
    e1 = zeros(len_x,1);
    e1(1) = 1;
    w = (x+niu*e1) * sqrt(k) / (x(1)+niu);
    sigma = (x(1)+niu) / (niu*k);
    delta = sigma * w * w';
end

function [delta] = householder2(x, len_x)
    x_norm = norm(x, 'fro');
    niu = sign(real(x(1))) * x_norm;
    e1 = zeros(len_x,1);
    e1(1) = 1;
    w = (x+niu*e1) / (x(1)+niu);
    sigma = (x(1)+niu) / niu;
    delta = sigma * w * w';
end