function [U,S,V]=jacobi_svd(A)
    % [U S V]=jacobi_svd(A)
    % A is original matrix
    % Singular values come back in S (diag matrix)
    % orig matrix = U*S*V’
    %
    % One-sided Jacobi algorithm for SVD

    TOL=1.e-8;
    MAX_STEPS=40;
    n=size(A,1);
    U=A;
    V=eye(n);
    for steps=1:MAX_STEPS
        converge=0;
        s = col_square_sum(U)
        for j=2:n
            for k=1:j-1
                alpha=s(k)
                beta=s(j)
                gamma=0
                for i=1:n
                    gamma=gamma+U(i,k)*U(i,j)
                end
                converge=max(converge,abs(gamma)/sqrt(alpha*beta));
                if gamma ~= 0
                    zeta=(beta-alpha)/(2*gamma);
                    t=sign(zeta)/(abs(zeta)+sqrt(1+zeta^2));
                else
                    t=0;
                end
                c=1/sqrt(1+t^2)
                s=c*t
                % update columns k and j of U
                t=U(:,k);
                U(:,k)=c*t-s*U(:,j);
                U(:,j)=s*t+c*U(:,j);
                % update matrix V of right singular vectors
                t=V(:,k);
                V(:,k)=c*t-s*V(:,j);
                V(:,j)=s*t+c*V(:,j);
            end
        end
        if converge < TOL
            break;
        end
    end

    if steps >= MAX_STEPS
        error(’jacobi_svd failed to converge!’);
    end
    % the singular values are the norms of the columns of U
    % the left singular vectors are the normalized columns of U
    for j=1:n
        singvals(j)=norm(U(:,j));
        U(:,j)=U(:,j)/singvals(j);
    end
    S=diag(singvals);

function s = col_square_sum(U)
    [n,m] = size(U);
    s = zeros(m,1);
    for i = 1:m
        sum = 0;
        for j = 1:n
            sum = sum + U(j,i)^2;
        end
        s(i) = sum;
    end


