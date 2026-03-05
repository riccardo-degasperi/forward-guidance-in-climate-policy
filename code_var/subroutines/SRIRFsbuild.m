function IRF = SRIRFsbuild(beta,B0,h,Ns,dims)

% Unpack
N  = dims.N;
p  = dims.p;
m  = dims.m;
mc = dims.mc;

% Containers
IRF = cell(N,Ns);
Yh  = NaN(h+1,N*p);

% Generate IRFs
for i = 1:Ns               %loop over shocks
Yh(1,:) = zeros(1,N*p);
Yh(1,1:N) = B0(:,i)';

    if h == 0
        
        for ii = 1:N
        IRF{ii,i} = Yh(:,ii);
        end
        
    else

        for j = 1:h       %loop over horizons
        Yh(j+1,:) = [Yh(j,:)*beta(1:end-m-mc,:) Yh(j,1:N*(p-1))];
        end

        for ii = 1:N
        IRF{ii,i} = Yh(:,ii);
        end
    
    end

end


