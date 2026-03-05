function [checkall,checkall_flip] = SRcheck(SR,IRF)

[N,~,Ns] = size(SR);
checkall      = ones(N,Ns);
checkall_flip = ones(N,Ns);

for j = 1:Ns
    for i = 1:N
        if SR(i,1,j) ~= 0
            if SR(i,3,j) == 1      %positive
                
                % Check signs
                check = IRF{i,j}(SR(i,1,j):SR(i,2,j),1) > 0;
                checkall(i,j) = min(check);
                
                % Check flipped signs
                check_flip = IRF{i,j}(SR(i,1,j):SR(i,2,j),1) < 0;
                checkall_flip(i,j) = min(check_flip);
                
            elseif SR(i,3,j) == -1 %negative
                
                % Check signs
                check = IRF{i,j}(SR(i,1,j):SR(i,2,j),1) < 0;
                checkall(i,j) = min(check);
                
                % Check flipped signs
                check_flip = IRF{i,j}(SR(i,1,j):SR(i,2,j),1) > 0;
                checkall_flip(i,j) = min(check_flip);
                
            end
        end
    end
end