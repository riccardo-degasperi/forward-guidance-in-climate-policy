function Q = gramschmidt(A)
% Gram-Schmidt Orthogonalization

    [m,n]=size(A);
    Q=zeros(m,n);
    R=zeros(n,n);
    for j=1:n
        v=A(:,j);
        for i=1:j-1
            R(i,j)=Q(:,i)'*v;
            v=v-R(i,j)*Q(:,i);
        end
        R(j,j)=norm(v);
        Q(:,j)=v/R(j,j);
    end
end

% ALTERNATIVE VERSION

% [~,n] = size(A);
% Q = zeros(n,n);
% tmp1 = zeros(n,n);   
% 
% for k = 1:n
% 
%     x = A(:,k);
% 
%     for i = 1:k-1
%         tmp1(:,i) = (Q(:,i)'*x)*Q(:,i);
%     end
%     tmp2 = x-sum(tmp1,2);
%     Q(:,k) = tmp2/sqrt(tmp2'*tmp2);
% end
