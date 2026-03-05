function out = BVARidentification(in,data,dims,modelOpt,dataOpt,miscOpt)
% Identifies a Bayesian VAR

% INPUT:
% - in       : output from BLPestimate.m
% - data     : output from makeXY;
% - dims     : output from makeXY;
% - modelOpt : structure containing model options;
% - dataOpt  : structure containing data options;
% - miscOpt  : structure containing misc options; 

% OUTPUT:
% - out      : structure containing estimation results and IRFs;


% r.degasperi@warwick.ac.uk
%-------------------------------------------------------------------------%

% Options
precision = 1.0000e-05;

%-------------------------------------------------------------------------%
% Unpacking
betas           = in.betas;
sigmas          = in.sigmas;
es              = in.es;
hpVAR           = in.hpVAR;
maxeig          = in.maxeig;
neig            = in.neig;
IRF             = in.IRF;
dates           = data.dates;
X               = data.X;
Y               = data.Y;
T               = dims.T;
N               = dims.N;
p               = dims.p;
m               = dims.m;
mc              = dims.mc;
h               = dims.h;
sim             = modelOpt.sim;
burnin          = modelOpt.burnin;
jump            = modelOpt.jump;
draws           = (sim-burnin)/jump;
identification  = modelOpt.identification;
varendo         = dataOpt.varendo;
shockVar        = dataOpt.shockVar;
shockSize       = dataOpt.shockSize;
shockSign       = dataOpt.shockSign;

% Load sign restrictions
if strcmp(modelOpt.identification,'jarocinski-karadi')
    SR          = modelOpt.SR;
    maxSRtries  = modelOpt.maxSRtries;
    if ~isempty(dataOpt.mrScript)
        MR      = modelOpt.MR;
        noMR    = 0;
    else
        noMR    = 1;
    end
end
if strcmp(modelOpt.identification,'iv')
    iv          = data.iv;
    ivDates     = data.ivDates;
end

% Check whether the shock is normalised or one SD
check1sd = strcmp('1sd',shockSize);


%-------------------------------------------------------------------------%
% Structural identification
disp('--> Recovering structural form...')

% Estimate impact matrix
switch identification
    
    case 'cholesky' %-----------------------------------------------------%
        
        B0nnr = NaN(N,N,draws);
        B0    = NaN(N,N,draws);
        
        for s = 1:draws
            
            % Cholesky decomposition
            B0nnr(:,:,s) = chol(nspd(sigmas(:,:,s)),'lower');

            % Unit normalisation
            B0(:,:,s) = B0nnr(:,:,s)*diag(1./diag(B0nnr(:,:,s)));
            
        end
        
        
    case 'iv' %-----------------------------------------------------------%
        
        % Initialise containers
        B0nnr = NaN(N,N,draws);
        
        % Match length of residual series to IV
        [es2,iv2,~,~,~,~] = ivMatch(p,dates,ivDates,es,iv);

        % Estimate identified columns of B0
        [B0_nr_part,B0_nnr_part,~,res] = ivB0(es2,iv2,sigmas,dims,dataOpt);

        % First-stage F-stat: median and bands
        Fdistr = [quantile(res.Fstat,0.05)' quantile(res.Fstat,0.5)' quantile(res.Fstat,0.95)'];
        Frdistr = [quantile(res.Frobust,0.05)' quantile(res.Frobust,0.5)' quantile(res.Frobust,0.95)'];

        % Correlation between (principal components of) proxies and shocks (see Mertens & Ravn, 2013, p.8)
        RMdistr = [quantile(res.corrIV,0.05)' quantile(res.corrIV,0.5)' quantile(res.corrIV,0.95)']; 

        % Get frequentist 1st stage F statistic
        if T-p > N*p+m+mc
            resid = Y-X*(X\Y);
            sigmaf = resid'*resid/(T-p-N*p-m-mc);
            [resid1s,iv1s,~,~,~,~] = ivMatch(p,dates,ivDates,resid,iv);
            [~,~,~,res_freq] = ivB0(resid1s,iv1s,sigmaf,dims,dataOpt);
            Fstat = res_freq.Fstat;
            Frobust = res_freq.Frobust;
        else
            Fstat = 'NA';
            Frobust = 'NA';
            warning('Frequentist F-stat cannot be computed because T-p < N*p+m')
        end
        
        % Obtain remaining columns of B0 (without structural meaning)
        % (Venditti & Veronese, 2021)

        % Specify number of shocks
        Ns = numel(shockVar);
        
        for s = 1:draws

            %1) Reorder policy variables first
            sigma = sigmas(:,:,s);
            iP = ismember(varendo,shockVar); %position of policy variables
            b1 = [B0_nnr_part(iP,iP,s);B0_nnr_part(~iP,iP,s)];
            Sig = [sigma(iP,iP) sigma(iP,~iP); sigma(~iP,iP) sigma(~iP,~iP)];
            
            %2) Compute Chol(Sigma) and omega1hat:
            Sigtr = chol(0.5*(Sig + Sig'))';
            Omega1hat = Sigtr\b1;
            % NOTE: We use the non-normalized matrix of impacts

            %3) Draw Omega (nxn):
            x = randn(N);
            Omega_draw = getqr(x);

            %4) Replace omega1hat in Omega:
            Omega_draw(:,1:Ns) = Omega1hat;

            %5) Gram-Schmidt orthogonalization
            Omegastar = gramschmidt(Omega_draw);
            % NOTE: if the additional restrictions used to distinguish the
            % shocks when you identify more than one shock with more than
            % one proxy are such that the identified shocks are not
            % orthogonal, this procedure will ortogonalise them,
            % potentially altering the impact coefficients on the shocks
            % ordered after the first one. If you use Choleski, the
            % identified columns will remain the same.
            if Ns>1 && sum(sum(Omega1hat(:,2:Ns)-Omegastar(:,2:Ns))) > precision
                error('ERROR: Gram-Schmidt is changing the proxy-identified columns.')
            end

            %6) Compute impact matrix B
            tmpB0nnr = Sigtr*Omegastar;

            %7) Check that the covariance matrix is maintained
            Sig_check = (Sigtr*Omegastar)*(Sigtr*Omegastar)';
            test = abs(Sig-Sig_check) > precision;
            if sum(sum(test))>0
                error('ERROR: The orthogonalized covariance differs from the input covariance.')
            end

            %8) Reorder columns
            B0nnr(iP,iP,s)   = tmpB0nnr(1:Ns,1:Ns);
            B0nnr(~iP,iP,s)  = tmpB0nnr(Ns+1:end,1:Ns);
            B0nnr(iP,~iP,s)  = tmpB0nnr(1:Ns,Ns+1:end);
            B0nnr(~iP,~iP,s) = tmpB0nnr(Ns+1:end,Ns+1:end);

        end
        
        % Store normalised partial impact matrix for later use
        B0 = B0_nr_part;


    case 'jarocinski-karadi' %--------------------------------------------%
        
        % Patch to compute normalised B0s even when shockSize = '1sd' 
        if strcmp('1sd',shockSize)
            shockSize = 1;
        end

        % Initialize containers for impact matrices
        B0nnr = zeros(N,N,draws);
        B0 = zeros(N,N,draws);
        QB01 = zeros(N,N);

        % Maximal length of IRF function to be checked
        nstepSR = max(max(SR(:,2,:)))-1;  %0:impact, 1:1 horizon ahead, ...
        if noMR
            nstepMR = 0;
        else
            nstepMR = max(max(MR(:,2,:)))-1;  %0:impact, 1:1 horizon ahead, ...
        end
        nsteps = max([nstepSR nstepMR]);  %0:impact, 1:1 horizon ahead, ...

        jj = 1; % accepted draws
        kk = 0; % total draws
        ww = 1; % index for printing on screen
        while jj <= draws

            kk = kk+1;
            if kk/jj > maxSRtries
                error(['ERROR: max number of iterations reached (accepted draws: ',num2str(jj),'). Try different sign/magnitude restrictions.'])
            end

            % Find position of shockVar
            iP = ismember(varendo,shockVar); %position of policy variables
            Ns = sum(iP);

            % Draw a random orthonormal matrix
            Q = eye(N);
            K = randn(Ns,Ns);
            Q0 = getqr(K);
            Q(1:Ns,1:Ns) = Q0;

            % Check precision
            if sum(sum(Q*Q'))> N + precision
                error('ERROR: Q*transpose(Q) is not equal to identity.')
            end

            % Draw beta and sigma from the posterior
            beta = betas(:,:,jj);
            sigma = sigmas(:,:,jj);

            % Reorder policy variables first
            Sig = [sigma(iP,iP) sigma(iP,~iP); sigma(~iP,iP) sigma(~iP,~iP)];

            % Compute candidate draw
            B0_ = chol(nspd(Sig));
            QB00 = B0_'*Q';

            % Reorder rows and columns of draw of impact matrix
            QB01(iP,iP)   = QB00(1:Ns,1:Ns);
            QB01(~iP,iP)  = QB00(Ns+1:end,1:Ns);
            QB01(iP,~iP)  = QB00(1:Ns,Ns+1:end);
            QB01(~iP,~iP) = QB00(Ns+1:end,Ns+1:end);

            % Reorder identified shocks first
            QB0_ = [QB01(:,iP) QB01(:,~iP)];
            
            z = 0;
            flag = 0;
            while z < 1
    
                % Compute IRFs only for the restricted periods
                IRF_draw = SRIRFsbuild(beta,QB0_,nsteps,Ns,dims);

                % Check if sign restrictions are satisfied
                [checkSR,checkSR_flip] = SRcheck(SR,IRF_draw);
                if noMR
                    checkMR = ones(size(checkSR));
                    checkMR_flip = ones(size(checkSR_flip));
                else
                    [checkMR,checkMR_flip] = MRcheck(MR,IRF_draw);
                end
                checkall = min(checkSR,checkMR);
                checkall_flip = min(checkSR_flip,checkMR_flip);
    
                % NOTE: Changing the sign of columns results in another 
                % orthonormal matrix (See Rubio-Ramirez, Waggoner & Zha, 2010)
                if min(min(checkall))==0

                    flip = sum(checkall,1)==0;
                    if sum(flip)>0 && flag==0
                        QB0_(:,flip)=-QB0_(:,flip);
                        flag = flag+1;
                    else
                        z = 1;
                    end
                else
                    z = 1;
                end
            end

            % If restrictions are satisfied store impact matrix
            if min(min(checkall))==1

                % Reorder columns of B0
                B0nnr(:,iP,jj)   = QB0_(:,1:Ns);
                B0nnr(:,~iP,jj)  = QB0_(:,Ns+1:end);
    
                % Normalize B0 relative to the first variable
                tmpb1 = QB0_(:,1:Ns)./QB0_(1,1:Ns)*shockSize;
                B0(:,iP,jj) = tmpb1(:,1:Ns);
                
                % Next iteration
                jj = jj+1;

                % Display number of iterations
                if jj == 100*ww
                    disp(['SR iteration: ',num2str(jj),'/',num2str(kk),' draws']);
                    ww = ww+1;
                end

            elseif min(min(checkall_flip))==1

                % Reorder columns of B0
                B0nnr(:,iP,jj)   = -QB0_(:,1:Ns);
                B0nnr(:,~iP,jj)  = -QB0_(:,Ns+1:end);
    
                % Normalize B0 relative to the first variable
                tmpb1 = QB0_(:,1:Ns)./QB0_(1,1:Ns)*shockSize;
                B0(:,iP,jj) = tmpb1(:,1:Ns);
                
                % Next iteration
                jj = jj+1;

                % Display number of loops
                if jj == 100*ww
                    disp(['SR iteration: ',num2str(jj),'/',num2str(kk),' draws']);
                    ww = ww+1;
                end
            end
        end


end


%-------------------------------------------------------------------------%
% Compute shock series
shock = zeros(T-p,numel(shockVar),draws);
iP = ismember(varendo,shockVar);
for s = 1:draws

    if check1sd

        % One SD normalisation
        shock(:,:,s) = (B0nnr(:,iP,s)'/sigmas(:,:,s)*es(:,:,s)')'; 
    
        % Alternative: one SD normalisation (as in Cesa-Bianchi)
        %shock_(:,:,s) = (B0nnr(:,:,s)\es(:,:,s)')';

    else

        % Unit normalisation
        shock(:,:,s) = (B0(:,iP,s)'/sigmas(:,:,s)*es(:,:,s)')'/(B0(:,iP,s)'/sigmas(:,:,s)*B0(:,iP,s));
    end
end


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
IRFs = IRFsbuild(betas,IRF,B0_,dims,shockSign);


%-------------------------------------------------------------------------%
% Get IRFs median, upper and lower bounds
IRFsulm = IRFbands(IRFs,N,h);


%-------------------------------------------------------------------------%
% Pack output
out.betas        = betas;
out.sigmas       = sigmas;
out.es           = es;
out.hpVAR        = hpVAR;
out.maxeig       = maxeig;
out.neig         = neig;
out.IRF          = IRF;
out.B0nnr        = B0nnr;
out.B0           = B0;
out.IRFs         = IRFs;
out.IRFsulm      = IRFsulm;
out.shock        = shock;
if strcmp(identification,'iv')
    out.Fstat    = Fstat;
    out.Frobust  = Frobust;
    out.Fdistr   = Fdistr;
    out.Frdistr  = Frdistr;
    out.RMdistr  = RMdistr;
end

end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% QR decomposition
function out = getqr(a)
    [q,r] = qr(a,0);
    out = q*diag(sign(diag(r)));

% Alternative version
% % (See Rubio-Ramirez, Waggoner & Zha, 2010)
% for i = 1:size(q,1)
%     if r(i,i) < 0
%         q(:,i) = -q(:,i);
%     end
% end
% out = q;

end
