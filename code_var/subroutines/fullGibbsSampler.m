function out = fullGibbsSampler(data,dims,hp,modelOpt,dataOpt,miscOpt)
% Full Gibbs Sampler Algorithm for Bayesian VARs

% INPUT:
% - data     : output from makeXY;
% - dims     : output from makeXY;
% - hp       : hyperparameters and related options;
% - modelOpt : structure containing model options;
% - dataOpt  : structure containing data options; 
% - miscOpt  : structure containing misc options;

% OUTPUT:
% - out      : structure containing estimation results and IRFs;

% riccardo.degasperi@bancaditalia.it
%-------------------------------------------------------------------------%

% Unpacking
X               = data.X;
Y               = data.Y;
endo            = data.endo;
Xex             = data.Xex;
T               = dims.T;
N               = dims.N;
p               = dims.p;
pexo            = dims.pexo;
m               = dims.m;
mc              = dims.mc;
h               = dims.h;
dropExplosive   = miscOpt.dropExplosive;
sim             = modelOpt.sim;
burnin          = modelOpt.burnin;
jump            = modelOpt.jump;
identification  = modelOpt.identification;
draws           = (sim-burnin)/jump;
GLP             = hp.GLP;
shockVar        = dataOpt.shockVar;
shockSize       = dataOpt.shockSize;
shockSign       = dataOpt.shockSign;

% Check whether the shock is normalised or one SD
check1sd = strcmp('1sd',shockSize);


%-------------------------------------------------------------------------%
% Estimate univariate AR(order) for each endogenous variable
hp.SS  = NaN(N,1);         %container for residual standard deviations
hp.SS2 = NaN(N,1);         %container for residual variances
order  = 1;                %we fit an AR(order) for the whole sample T

for i = 1:N

    % Build X and Y matrices
    X2  = [ones(T-order,1) lagX(endo(:,i),order)]; % (T-p)x(p+1) matrix (with constant)
    Y2  = endo(order+1:end,i);                     % (T-p)x1 matrix
    
    % Intermediate step
    XX2 = X2'*X2;
    XY2 = X2'*Y2;

    % Compute beta and std of residuals
    b2         = XX2\XY2;
    hp.e2(:,i) = Y2-X2*b2;
    hp.SS(i)   = std(Y2-X2*b2);
    hp.SS2(i)  = hp.SS(i).^2;
end

% Compute mean across first p observations for endogenous
Ybar = mean(endo(1:p,:),1)';

% Compute mean across first p observations for exogenous
if ~isempty(Xex)
    Xbar = mean(Xex(1:p-pexo,:),1);
else
    Xbar = [];
end

%-------------------------------------------------------------------------%
% Optimal priors
hpVAR = hp;   %initialise hpVAR

if GLP
%disp('--> Computing optimal priors...')
    
    [hpglp,csmin,postMode] = GLPoptimalPriors(Y,X,Ybar,Xbar,hp,dims);

    % Update hyperparameters
    f     = fieldnames(hpglp);
    for i = 1:length(f)
        hpVAR.(f{i}) = hpglp.(f{i});
    end

end

% Generate dummy observations
[Yd,Xd] = dummyObs(Ybar,Xbar,hpVAR,dims);

% Combined X and Y matrices
Yd = [Y; Yd];
Xd = [X; Xd];

% Estimate reduced-form VAR
beta0    = Xd\Yd;
e0      = Yd-Xd*beta0;            %including dummies
[TT,NN] = size(Xd);
d       = N+2;
sigma0   = (e0'*e0)/(TT-NN+d+1);  %same Dof correction as GLP

%---------------------------------------------------------------------%
% Sampler

% Initialise containers
betas     = nan(N*p+m+mc,N,draws);
sigmas    = nan(N,N,draws);
es        = nan(T-p,N,draws);
B0nnr     = nan(N,N,draws);
B0        = nan(N,N,draws);
shock     = nan(T-p,numel(shockVar),draws);
maxeig    = nan(sim,1);
neig      = 0;

k = 1;                                          %initialise index for storage matrices
s = 1;                                          %initialise count
w = 1;                                          %initialise count

CinvS  = chol(inv((TT-NN+d+1)*nspd(sigma0)));   %removing DoF correction from sigma
CinvXX = chol(inv(Xd'*Xd));

while k <= draws

    % Display advancement
    if k == 10*w
        if exist('txt','var')
            fprintf(repmat('\b', 1, numel(txt)));
        end
        disp(' ')
        txt = ['    Accepted draws: ',num2str(k),'/',num2str(s),' (of ',num2str(draws),' required) [',datestr(datetime),']\n'];
        fprintf(txt);
        w = w+1;
    end

    %Get Psi
    Xg = randn(TT+2-(N*p+m+mc),N)*CinvS;           % vec(Xg) ~ N(vec(0),kron(inv(TT*sigma),I))
    sigma = inv(Xg'*Xg);                        % (Xg'Xg) ~ W(inv(TT*sigma),TT+2-(N*p+1))
    % NOTE: Dof same as Banbura, Giannone & Reichlin (2010).
    %       TT+2-(N*p+m) = T+N+2 (= df in GLP)
    
    %Get Cholesky factor to use in MN below
    Csigmas = chol(sigma);

    %Draw from posterior of beta (matrixvariate)
    tmp  = randn(size(beta0));
    beta = beta0 + CinvXX'*tmp*Csigmas;         % vec(betag) ~ N(vec(beta),kron(Sg,invXX))
    
    % Check maximum eigenvalue
    A = NaN(N*p);                               %create companion matrix
    A(1:N,:) = beta(1:end-m-mc,:)';
    A(N+1:end,:) = [eye(N*(p-1)) zeros(N*(p-1),N)];
    maxeig(s,1) = max(abs(eig(A)));             %compute maximum eigenvalue
    if maxeig(s,1) > 1
        neig = neig+1;
    end

    % Get residuals
    eps = Y - X*beta;

    % Burnin phase
    if s <= burnin
        s = s+1;
    
    % If we keep explosive draws
    elseif ~dropExplosive && s > burnin
        
        if mod(s,jump) == 0

            % Estimate impact matrix
            in.betas   = beta;
            in.sigmas  = sigma;
            in.es      = eps;
            [IDX,flag] = fullGibbsIdentification(in,data,dims,modelOpt,dataOpt);

            % Skip draw if the covariance from gram-schmidt fails
            if flag == 0

                B0nnr(:,:,k)     = IDX.B0nnr;
                B0(:,:,k)        = IDX.B0;
                shock(:,:,k)     = IDX.shock;
                betas(:,:,k)     = beta;               %VAR coefficients
                sigmas(:,:,k)    = sigma; 
                es(:,:,k)        = Y - X*beta;         %residuals for structural identification

                k = k+1;                               %update index for storage matrices
            end
        end
        
        s = s+1;

    % If we discard explosive draws
    elseif dropExplosive && s > burnin && maxeig(s,1) < 1
    
        if mod(s,jump) == 0

            % Estimate impact matrix
            in.betas   = beta;
            in.sigmas  = sigma;
            in.es      = eps;
            [IDX,flag] = fullGibbsIdentification(in,data,dims,modelOpt,dataOpt);
    
            % Skip draw if the covariance from gram-schmidt fails
            if flag == 0

                B0nnr(:,:,k)     = IDX.B0nnr;
                B0(:,:,k)        = IDX.B0;
                shock(:,:,k)     = IDX.shock;
                betas(:,:,k)     = beta;               %VAR coefficients
                sigmas(:,:,k)    = sigma; 
                es(:,:,k)        = Y - X*beta;         %residuals for structural identification

                k = k+1;                               %update index for storage matrices
            end
        end
        
        s = s+1;
        
    end
end

% Display number of explosive draws
if ~dropExplosive
disp(['    Number of explosive draws: ' num2str(neig)])
end
disp(' ')


%-------------------------------------------------------------------------%
% Generate Structural-form Impulse Response Functions
disp('--> Generating structural-form IRFs...')

% If shockSize = '1sd', use non-normalised impact matrix to generate IRFs
if check1sd
    B0_ = B0nnr;
else
    B0_ = B0;
end

% Generate structural impulse responses
IRF = []; %unused
IRFs = IRFsbuild(betas,IRF,B0_,dims,shockSign);



%-------------------------------------------------------------------------%
% Get IRFs median, upper and lower bounds
IRFsulm = IRFbands(IRFs,N,h);


%-------------------------------------------------------------------------%
% Pack output
out.betas         = betas;
out.sigmas        = sigmas;
out.es            = es;
out.hpVAR         = hpVAR;
out.maxeig        = maxeig;
out.neig          = neig;
out.B0nnr         = B0nnr;
out.B0            = B0;
out.shock         = shock;
out.IRFs          = IRFs;
out.IRFsulm       = IRFsulm;

if strcmp(identification,'iv')
    out.Fstat     = IDX.Fstat;
    out.Frobust   = IDX.Frobust;
end

