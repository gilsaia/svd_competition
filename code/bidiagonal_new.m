function [U,B,V] = bidiagonal_new(A,r)
    [m,n]=size(A);
    B=A;
    U=eye(m);
    V=eye(n);
    for j=1:n
        [sigma,w] = householderq(B(j:end,j),m-j+1);  

        wt=sigma*w;
        bt=wt'*B(j:end,j:end);      
        B(j:end,j:end) = B(j:end,j:end) - w*bt;
        % wt=scalemat(sigma,w);
        % bt=matmul(wt',B(j:end,j:end));
        % B(j:end,j:end)=B(j:end,j:end)-vecmulvectomat(w,bt);

        ut=U(:,j:end)*w;
        U(:,j:end)=U(:,j:end)-ut*wt';
        % ut=matmul(U(:,j:end),w);
        % U(:,j:end)=U(:,j:end)-vecmulvectomat(ut,wt');

        if j <= n-1
            [sigma,w] = householderp(B(j,j+1:end)',n-j);

            wt=sigma*w;
            bt=B(j:end,j+1:end)*wt;
            B(j:end,j+1:end)=B(j:end,j+1:end)-bt*w';
            % wt=scalemat(sigma,w);
            % bt=matmul(B(j:end,j+1:end),wt);
            % B(j:end,j+1:end)=B(j:end,j+1:end)-vecmulvectomat(bt,w');

            vt=w'*V(j+1:end,:);
            V(j+1:end,:)=V(j+1:end,:)-wt*vt;
            % vt=matmul(w',V(j+1:end,:));
            % V(j+1:end,:)=V(j+1:end,:)-vecmulvectomat(wt,vt);
        end
    end

end

function [sigma,w] = householderp(x,len_x)
    x_norm = norm(x, 'fro');
    niu = sign(real(x(1))) * x_norm;
    k = (abs(real(x(1)))+x_norm) / x_norm;
    xt=x;
    xt(1)+=niu;
    w = xt * sqrt(k) / (x(1)+niu);
    sigma = (x(1)+niu) / (niu*k);
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

function [sigma,w] = householderq(x,len_x)
    x_norm = norm(x, 'fro');
    niu = sign(real(x(1))) * x_norm;
    xt=x;
    xt(1)+=niu;
    w = xt / (x(1)+niu);
    sigma = (x(1)+niu) / niu;
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