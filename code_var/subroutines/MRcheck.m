function [checkall,checkall_flip] = MRcheck(MR,IRF)

[N,~,Ns] = size(MR);
checkall      = ones(N,Ns);
checkall_flip = ones(N,Ns);

for j = 1:Ns
    for i = 1:N
        if MR(i,1,j) ~= 0
            if MR(i,4,j) == 1      %positive

                if MR(i,5,j) == i
                
                % Check signs
                check = IRF{i,j}(MR(i,1,j):MR(i,2,j),1) > MR(i,3,j);
                checkall(i,j) = min(check);
                
                % Check flipped signs
                check_flip = IRF{i,j}(MR(i,1,j):MR(i,2,j),1) < -MR(i,3,j);
                checkall_flip(i,j) = min(check_flip);

                else

                pos = MR(i,5,j);

                % Check signs
                check = IRF{i,j}(MR(i,1,j):MR(i,2,j),1) > MR(i,3,j).*IRF{pos,j}(MR(i,1,j):MR(i,2,j),1);
                checkall(i,j) = min(check);
                
                % Check flipped signs
                check_flip = IRF{i,j}(MR(i,1,j):MR(i,2,j),1) < -MR(i,3,j).*IRF{pos,j}(MR(i,1,j):MR(i,2,j),1);
                checkall_flip(i,j) = min(check_flip);

                end
                
            elseif MR(i,4,j) == -1 %negative

                if MR(i,5,j) == i
                
                % Check signs
                check = IRF{i,j}(MR(i,1,j):MR(i,2,j),1) < MR(i,3,j);
                checkall(i,j) = min(check);
                
                % Check flipped signs
                check_flip = IRF{i,j}(MR(i,1,j):MR(i,2,j),1) > -MR(i,3,j);
                checkall_flip(i,j) = min(check_flip);

                else

                pos = MR(i,5,j);

                % Check signs
                check = IRF{i,j}(MR(i,1,j):MR(i,2,j),1) < MR(i,3,j).*IRF{pos,j}(MR(i,1,j):MR(i,2,j),1);
                checkall(i,j) = min(check);
                
                % Check flipped signs
                check_flip = IRF{i,j}(MR(i,1,j):MR(i,2,j),1) > -MR(i,3,j).*IRF{pos,j}(MR(i,1,j):MR(i,2,j),1);
                checkall_flip(i,j) = min(check_flip);

                end
                
            end
        end
    end
end