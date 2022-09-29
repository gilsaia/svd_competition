function [U,S,V]=jacobi_svd(A)
    % [U S V]=jacobi_svd(A)
    % A is original matrix
    % Singular values come back in S (diag matrix)
    % orig matrix = U*S*V’
    % One-sided Jacobi algorithm for SVD

    % !!todo: check the dim of U

    TOL=1.e-8;
    MAX_STEPS=40;
    n=size(A,2);
    U=A;
    V=eye(n);
    for steps=1:MAX_STEPS
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
                    U(i,k)=c*U(i,k)-st*U(i,j);
                    U(i,j)=st*U(i,k)+c*U(i,j);
                end
                % update matrix V of right singular vectors
                for i=1:n
                    V(i,k)=c*V(i,k)-st*V(i,j);
                    V(i,j)=st*V(i,k)+c*V(i,j);
                end
            end
        end
        % the singular values are the norms of the columns of U
        % the left singular vectors are the normalized columns of U
        UT=U;
        for j=1:n
            singvals(j)=norm(U(:,j),'fro');
            for i=1:n
                UT(i,j)=UT(i,j)/singvals(j);
            end
        end
        S=diag(singvals);
        AT=UT*S*V;
        Remain=A-AT;
        disp(sprintf('Step:%d Coverge:%f,Norm:%f',steps,converge,norm(Remain,'fro')/norm(A,'fro')));
        if converge < TOL
            break;
        end
    end

    if steps >= MAX_STEPS
        error('jacobi_svd failed to converge!');
    end
    % the singular values are the norms of the columns of U
    % the left singular vectors are the normalized columns of U
    for j=1:n
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


