function [U,S,V]=jacobi(A,V,n,r)
    % [U S V]=jacobi_svd(A)
    % A is bidiagonal matrix
    % Singular values come back in S (diag matrix)
    % orig matrix = U*S*V'
    % One-sided Jacobi algorithm for SVD


    TOL=1e-8;
    MAX_STEPS=40;
    U=A;
    for steps=1:MAX_STEPS
        tic;
        converge=0;
        s = col_square_sum(U, n);
        for j=2:n
            for k=1:j-1
                alpha=s(k);
                beta=s(j);
                gamma=0;
                for i=1:n
                    gamma=gamma+U(i,k)*U(i,j);
                end
                converge=max(converge,abs(gamma)/sqrt(alpha*beta));
                if gamma ~= 0
                    zeta=(beta-alpha)/(2*gamma);
                    t=sign(zeta)/(abs(zeta)+sqrt(1+zeta^2));
                else
                    t=0;
                end
                c=1/sqrt(1+t^2);
                st=c*t;
                % update columns k and j of U
                for i=1:n
                    tmp=U(i,k);
                    U(i,k)=c*U(i,k)-st*U(i,j);
                    U(i,j)=st*tmp+c*U(i,j);
                end
                % update matrix V of right singular vectors
                for i=1:n
                    tmp=V(i,k);
                    V(i,k)=c*V(i,k)-st*V(i,j);
                    V(i,j)=st*tmp+c*V(i,j);
                end
            end
        end
        steptime=toc;
        disp(sprintf('Step:%d Coverge:%f Time:%f',steps,converge,steptime));
        if converge < TOL
            break;
        end
    end

    if steps >= MAX_STEPS
        error('jacobi_svd failed to converge!');
    end
    % the singular values are the norms of the columns of U
    % the left singular vectors are the normalized columns of U
    singvals=zeros(n,1);
    for j=1:r
        singvals(j)=norm(U(:,j),'fro');
        for i=1:n
            U(i,j)=U(i,j)/singvals(j);
        end
    end
    S=diag(singvals);

function s = col_square_sum(U, n)
    s = zeros(n,1);
    for i = 1:n
        s(i) = norm(U(:,i), 'fro')^2;
    end


