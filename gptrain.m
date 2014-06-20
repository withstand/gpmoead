function gp = gptrain( data, y, oldgp )
%GP_HYP_EST Summary of this function goes here
%   Detailed explanation goes here
gp = oldgp;
[~, dimension] = size(data);

Delta = getfieldwithdefault(oldgp, 'Delta', 0);

lhp = getfieldwithdefault(oldgp,'lhp',[-6*ones(1,dimension) ones(1,dimension)]);
uhp = getfieldwithdefault(oldgp,'uhp',[6*ones(1,dimension) 2*ones(1,dimension)]);

hp0 = getfieldwithdefault(oldgp,'hp0',...
    [log10(getfieldwithdefault(oldgp,'theta', ones(1,dimension))), ...
    getfieldwithdefault(oldgp,'p', 2*ones(1,dimension))]);
%     [ ]);

if sum(Delta)~= 0  %noisy
    %     2 * dimension + 1
    s2 = std(y) * std(y);
    lsigma2 =  getfieldwithdefault(oldgp,'lsigma2', 1e-5*s2);
    usigma2 = getfieldwithdefault(oldgp,'lsigma2', 1e5*s2);
    sigma20 = getfieldwithdefault(oldgp,'lsigma20', ...
        getfieldwithdefault(oldgp,'sigma2',s2));
    lhp = [lsigma2 lhp];
    uhp = [usigma2 uhp];
    hp0 = [sigma20 hp0];
else  %noisy-free
    %     2 * dimension
end

refresh = getfieldwithdefault(oldgp,'refresh',0);
if refresh
    strDisplay = 'iter';
else
    strDisplay = 'none';
end
algorithm = getfieldwithdefault(oldgp,'algorithm', 'pd');
if strcmp(algorithm,'de')
    infPopsize = getfieldwithdefault(oldgp,'infPopsize', ...
        max(80,10*numel(lhp)));
    infGeneration = getfieldwithdefault(oldgp, 'infGeneration', ...
        max(200,50*numel(lhp)));
    [hp,min_nlml]  = optimizer_de(infPopsize, infGeneration);
elseif strcmp(algorithm,'pd')
    [hp,min_nlml] = fmincon(@nlml_pd, hp0,[],[],[],[], lhp,uhp, ...
        [],optimset('Display',strDisplay,...%'Algorithm','interior-point', ...
        'GradObj','on'));   
else
    error('Unrecognized algorithm %s', algorithm);
end

gp = hyper2opts(gp, hp);
gp.nlml = min_nlml;
% return hyperparameters as gpOpts struct
gp = gp_predict(data, y, gp);
gp.gp_func = @gpevaluate;

    function [hp,min_nlml] = optimizer_de(NP,itermax)
        %% Setup Optimization Algorithm to solve minimize minus-log-marginal-likelihood
        S_struct.F_VTR        = -inf;
        S_struct.FVr_minbound = lhp;
        S_struct.FVr_maxbound = uhp;
        S_struct.I_D          = numel(lhp);
        S_struct.I_NP         = NP; %max(80,10*hyperparameterSize);
        S_struct.I_itermax    = itermax; %;
        
        
        S_struct.F_weight     = 0.85;
        S_struct.F_CR         = 1;
        S_struct.I_strategy   = 2;
        
        S_struct.I_bnd_constr = 1;
        
        S_struct.I_refresh    = refresh;%getfieldwithdefault(gp_opts,'refresh',0);
        S_struct.I_plotting   = 0;
        
        [hp,s_v,~] = deopt(@nlmlde,S_struct);
        min_nlml = s_v.FVr_oa;
        
        
    end


%%
    function Obj = nlmlde(parameters,~)
        Obj.FVr_oa = nlml(parameters);
        Obj.I_nc = 0;
        Obj.FVr_ca = 0;
        Obj.I_no = 1;
    end

%%
    function likelihood = nlml(parameters)
        gpOpts_nlml = hyper2opts([], parameters);
        gp_nlml = gp_predict(data,y, gpOpts_nlml);
        likelihood = gp_nlml.nlml;
    end
%%
    function [likelihood,likelihood_pd] = nlml_pd(parameters)
        gpOpts_nlml = hyper2opts([], parameters);
        [gp_nlml,likelihood_pd] = gp_predict(data,y, gpOpts_nlml);
        likelihood = gp_nlml.nlml;
    end

%%
    function opts = hyper2opts(opts_orig, hyper)
        opts = opts_orig;
        opts.theta = 10.^hyper(end-2*dimension+1:end-dimension);
        opts.p = hyper(end-dimension+1:end);
        if sum(Delta)~=0
            opts.sigma2 = hyper(1);
        end
    end

    function yStar = gpevaluate(xStar)
        %% make prediction at specific point xStar
        
        [sampleSize, dimension] = size(gp.data);
        one = ones(sampleSize,1);
        [r,c] = size(xStar);
        if c == dimension
            newSites = r;
            xS =  xStar;
        elseif r ==dimension
            xS = xStar';
            newSites = c;
        end
        yStar = zeros(newSites,2);
        %   xS = reshape(xStar, 1, numel(xStar));
        Delta = getfieldwithdefault(gp,'Delta',0);
        if sum(Delta)~=0 %noisy
            k = gp.sigma2 * gp.corr(gp.data, xS, gp);
            for i=1:newSites
                invK_k = gp.L' \(gp.L \ k(:,i));
                hat_y = gp.mu + k(:,i)' * gp.invK_y_mu;
                hat_std2 = gp.sigma2 - k(:,i)' * invK_k +...
                    (1-one'* invK_k)^2 / (one'*gp.invK_one);
                if hat_std2 < 0
                    hat_std2 = 0;
                end
                hat_std = sqrt(hat_std2);
                yStar(i,:) = [hat_y, hat_std];
            end
        else
            k = gp.corr(gp.data, xS, gp);
            for i=1:newSites
                invK_k = gp.L' \(gp.L \ k(:,i));
                hat_y = gp.mu + k(:,i)' * gp.invK_y_mu;
                hat_std2 = gp.sigma2 * (1 - k(:,i)' * invK_k +...
                    (1-one'* invK_k)^2 / (one'*gp.invK_one));
                if hat_std2 < 0
                    hat_std2 = 0;
                end
                % Return std derivation, not the square of it
                hat_std = sqrt(hat_std2);
                yStar(i,:) = [hat_y, hat_std];
            end
        end
    end

    function [gp,glml_pd] = gp_predict(data,y, oldGp)
        
        gp = oldGp;
        % helper function to get default correlation function
        if nargin == 0
            gp.corr = @expCorr;
            return;
        end
%         narginchk(3,3);
        
        % get default correlation function, or use provided correlation function
        corrFunc = getfieldwithdefault(oldGp,'corr',@expCorr);
        
        [sampleSize, dimension] = size(data);
        assert(sampleSize==numel(y));
        
        % Delta must be provided or the problem will be treated as noisy free
        Delta = getfieldwithdefault(oldGp,'Delta',0);
        if numel(Delta)==1
            Delta = Delta * ones(sampleSize,1);
        end
        
        
        
        %% trainning
        R = corrFunc(data, data, oldGp);
        
        if sum(Delta)~=0   % noisy
            assert(isfield(oldGp, 'sigma2'));
            sigma2_est = oldGp.sigma2;
            K = sigma2_est * R + diag(Delta);
        else % noisy-free
            K = R;
        end
        
        one = ones(sampleSize,1);
        [L,q] = chol(K,'lower');
        gp.modify_K = (q>0);
        gp.modify = 0;
        while (q>0)
            % warning('FCMCGP:covarianceMatrixNonPSD', 'K are not positive definite.')
            K = K + eye(sampleSize)*eps*bitshift(1,gp.modify);
            gp.modify = gp.modify + 1;
            [L,q] = chol(K,'lower');
        end
        
        
        invK_y = L' \ (L \ y);
        invK_one = L' \ (L \ one);
        mu_est = one' * invK_y / (one' * invK_one);
        y_mu = y - mu_est;
        detK = prod(diag(L));
        % test of zero det
        gp.zero_det_K =  (detK < realmin);
        invK_y_mu = L' \ (L \ y_mu);
        y_mu_invK_y_mu = y_mu' * invK_y_mu;
        if sum(Delta)~=0   % noisy
            negativeLML = y_mu_invK_y_mu ...
                + sampleSize * log(2*pi);
        else % noisy-free
            sigma2_est =  y_mu_invK_y_mu / sampleSize;
            negativeLML = y_mu_invK_y_mu / sigma2_est ...
                + sampleSize * log(sigma2_est*2*pi);
        end
        % if detK == 0
        % negativeLML --> -Inf
        % which is not acceptable in optimization
        % try to solve this by
        if ~gp.zero_det_K
            negativeLML = negativeLML + log(detK);
        end
        
        negativeLML = 0.5 * negativeLML;
        % noisy free get an estimation of sigma
        % if sum(Delta)==0 || 1
        % end
        
        
        if nargout >= 2
            if sum(Delta)~=0  % must be noisy, or it doesnot make sense
                glml_pd = zeros(2*dimension+1,1);
                
                pdK_sigma2 = R;
                glml_pd(1) = 0.5 * trace(L'\(L\pdK_sigma2)) - ...
                    0.5 * y_mu' * (L'\ (L\ (pdK_sigma2 * (L'\ (L\ y_mu)))));
                start = 1;
            else
                glml_pd = zeros(2*dimension,1);
                start = 0;
            end
            for k = 1:dimension
                xk = data(:,k);
                d_k = reshape(abs(xk(repmat((1:sampleSize)', sampleSize,1))- ...
                    xk(reshape(repmat(1:sampleSize, sampleSize,1), sampleSize*sampleSize,1))),...
                    sampleSize, sampleSize) / oldGp.theta(k);
                Kinv = L'\(L\eye(sampleSize));
                
                if sum(Delta)~=0
                    pdK_theta_k = sigma2_est *  oldGp.p(k) / oldGp.theta(k) * ...
                        (d_k .^ oldGp.p(k)) .* R;
                    pdK_p_k = - sigma2_est *  oldGp.p(k) / oldGp.theta(k) * ...
                        log(d_k) .* (d_k .^ oldGp.p(k)) .* R;
                    %                     pdK_p_k(((1:sampleSize)'-1) * sampleSize + (1:sampleSize)') = 0;
                    pdK_p_k(d_k==0)=0;
                    glml_pd(start+k) = 0.5 * trace(Kinv * pdK_theta_k) - ...
                        0.5 * y_mu' * (Kinv * (pdK_theta_k * (Kinv * y_mu)));
                    glml_pd(start+k+dimension) = 0.5 * trace(Kinv * pdK_p_k) - ...
                        0.5 * y_mu' * (Kinv * (pdK_p_k * (Kinv * y_mu)));
                else
                    pdK_theta_k = oldGp.p(k) / oldGp.theta(k) * ...
                        (d_k .^ oldGp.p(k)) .* R;
                    pdK_p_k = - oldGp.p(k) / oldGp.theta(k) * ...
                        log(d_k) .* (d_k .^ oldGp.p(k)) .* R;
                    %                   pdK_p_k(((1:sampleSize)'-1) * sampleSize + (1:sampleSize)') = 0;
                    pdK_p_k(d_k==0)=0;
                    
                    glml_pd(start+k) = 0.5 * trace(Kinv *pdK_theta_k) - ...
                        0.5/ sigma2_est * y_mu' * (Kinv * (pdK_theta_k * (Kinv * y_mu)));
                    glml_pd(start+k+dimension) = 0.5 * trace(Kinv *pdK_p_k) - ...
                        0.5 / sigma2_est * y_mu' * (Kinv * (pdK_p_k * (Kinv * y_mu)));
                end
            end
            if sum(isinf(glml_pd))~=0 || sum(isnan(glml_pd))~=0
                disp('eeori');
            end
        end
        
        gp.data = data;
        gp.y = y;
        gp.K = K;
        gp.L = L;
        gp.invK_y = invK_y;
        gp.invK_one = invK_one;
        gp.mu = mu_est;
        gp.sigma2 = sigma2_est;
        gp.detK = detK;
        gp.invK_y_mu = invK_y_mu;
        gp.y_mu_invK_y_mu = y_mu_invK_y_mu;

        gp.nlml = negativeLML;        
        gp.corr = corrFunc;
%         gp.gpOpts = gpOpts;
        
        
        
        
        
        function r = expCorr(x1,x2, opts)
            [s1, d1] = size(x1);
            [s2, d2] = size(x2);
            assert(d1==d2);
            theta = reshape(getfieldwithdefault(opts,'theta',ones(1,d1)),1,d1);
            p = reshape(getfieldwithdefault(opts,'p',2 * ones(1,d1)),1,d1);
            
            assert(d1 == numel(theta) && d1 == numel(p));
            s1s2 = s1*s2;
            r = reshape(exp(-sum((abs(x1(repmat((1:s1)', s2,1),:)- ...
                x2(reshape(repmat(1:s2, s1,1), s1s2,1),:))./ ...
                theta(ones(s1s2,1),:)) .^ p(ones(s1s2,1),:),2)),s1,s2);
            

            
        end
    end

end





