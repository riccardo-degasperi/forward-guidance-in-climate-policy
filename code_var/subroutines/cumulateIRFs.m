function out = cumulateIRFs(IRFs,iFD,dims)
%cumulateIRFs cumulates structural IRFs for differenced variables

% Unpack
N = dims.N;
h = dims.h;


for ii = 1:N      %shocks
    for i = 1:N   %variables

        if iFD(i) == 1

        tmp = IRFs{i,ii};

        IRFc = tmp;
        for j = 2:h+1
            IRFc(j,:) = IRFc(j,:) + IRFc(j-1,:);
        end

        IRFs{i,ii} = IRFc;

        end

    end
end

out = IRFs;

