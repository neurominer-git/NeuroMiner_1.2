function [sY, IN] = nk_PartialCorrelationsObj( Y, IN )
% =========================================================================
% FORMAT [sY, IN] = nk_PartialCorrelationsObj( Y, IN )
% =========================================================================
% Remove nuisance effects IN.G from Y (opt. using predefined estimators)
%
% I\O Arguments:
% -------------------------------------------------------------------------
% Y                 : M cases x N features data matrix
% IN.G              : The covariate(s) to regress out from Y
% IN.nointercept    : Include an intercept in the model or not
% IN.subgroup       : Index vector of cases in Y to compute the beta(s) from
% IN.beta           : The estimated beta coefficients
% IN.revertflag     : Increase (=1) or attenuate (=0) IN.G effects
% IN.METHOD         : Partial correlations (=1), ComBat (=2), fastICA (=3)
% IN.featind        : Feature subspace to apply the correction to           
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 06 / 2023

method = @PartialCorrelationsObj; methodsel = 1;
if exist('IN','var') && ~isempty(IN) && isfield(IN,'METHOD') && IN.METHOD == 2
    method = @CombatObj; methodsel = 2;
elseif exist('IN','var') && ~isempty(IN) && isfield(IN,'METHOD') && IN.METHOD == 3
    method = @FastICAObj; methodsel = 3;
end
if isfield(IN,'copy_ts') && IN.copy_ts
    str = 2; off = 1;
else
    str = 1; off = 0;
end
% =========================== WRAPPER FUNCTION ============================ 
if iscell(Y) && exist('IN','var') && ~isempty(IN)
    sY = cell(1,numel(Y));
    for i = str:numel(Y) + off
        switch methodsel 
            case 1
                if isfield(IN,'beta') && ~isempty(IN.beta)
                    IN.G = IN.TsCovars{i};
                else
                    IN.G = IN.TrCovars{i};
                end
            case 2
                if isfield(IN,'estimators') && ~isempty(IN.estimators)
                    IN.G = IN.TsCovars{i};
                    IN.M = IN.TsMod{i};
                else
                    IN.G = IN.TrCovars{i};
                    IN.M = IN.TrMod{i};
                end
            case 3
                if isfield(IN,'sigICs') && ~isempty(IN.sigICs)
                        IN.G = IN.TsCovars{i};
                else
                    if size(IN.TrCovars,2) > 1
                        IN.G = IN.TrCovars{i};
                    else
                        IN.G = IN.TrCovars; 
                    end
                end
        end
        [sY{i-off}, IN] = method (Y{i-off}, IN); 
    end
else
    switch methodsel 
        case 1
            if isfield(IN,'beta') && ~isempty(IN.beta)
                if iscell(IN.TsCovars)
                    IN.G = IN.TrCovars;
                else
                    IN.G = IN.TsCovars;
                end
            else
                IN.G = IN.TrCovars;
            end

        case 2
            if isfield(IN,'estimators') && ~isempty(IN.estimators)
                if iscell(IN.TsMod)
                    IN.G = IN.TrCovars;
                    IN.M = IN.TrMod;
                else
                    IN.G = IN.TsCovars;
                    IN.M = IN.TsMod;
                end
            else
                IN.G = IN.TrCovars;
                IN.M = IN.TrMod;
            end
        case 3
            if isfield(IN, 'sigICs') && ~isempty(IN.sigICs)
                if iscell(IN.TsCovars)
                    IN.G = IN.TrCovars;
                else
                    IN.G = IN.TsCovars;
                end
            else
                IN.G = IN.TrCovars;
            end
    end

    [ sY, IN ] = method( Y, IN );
end

% =========================================================================
function [Y, IN] = PartialCorrelationsObj( Y, IN )

if isempty(IN),eIN=true; else, eIN=false; end

if eIN|| ~isfield(IN,'G') || isempty(IN.G) 
    error('No covariates defined in parameter structure'), 
end

if eIN || (~isfield(IN,'nointercept') || isempty(IN.nointercept) || ~IN.nointercept ) 
     %Create intercept vecotr
    interceptflag = true;
    intercept = ones(size(IN.G,1),1);
    %Check if intercept is already included in matrix to avoid double
    %intercept removal
    if isfield(IN,'beta') && ~isempty(IN.beta)
        if size(IN.beta,1) == size(IN.G,2), interceptflag = false; end
    end
else
    interceptflag = false;
end

if interceptflag
    %fprintf(' ... adding intercept to covariate matrix')
    IN.G = [intercept IN.G];
end

if eIN || ~isfield(IN,'beta') || isempty(IN.beta) 
    if ~isfield(IN,'subgroup') || isempty(IN.subgroup)
        % Compute IN.beta from entire population
        IN.beta = pinv(IN.G) * Y; 
    else
        % Compute IN.beta from a subgroup of observations
        idxSubgroup = logical(IN.subgroup);
        IN.beta = pinv(IN.G(idxSubgroup,:)) * Y(idxSubgroup,:);
    end
end

if isfield(IN,'featind')
    featind = IN.featind;
else
    featind = true(1,size(Y,2));
end

if eIN || ~isfield(IN,'revertflag') || isempty(IN.revertflag) || ~IN.revertflag
    Y(:,featind) = Y(:,featind) - IN.G * IN.beta(:,featind);
else
    Y(:,featind) = Y(:,featind) + IN.G * IN.beta(:,featind);
    
end

% =========================================================================
function [ Y, IN ] = CombatObj( Y, IN )

if isempty(IN),eIN=true; else, eIN=false; end

if eIN|| ~isfield(IN,'G') || isempty(IN.G), error('No covariates defined in parameter structure'), end

% Combat needs data, batch variables and covariates transposed 
Y=Y'; G = IN.G'; if islogical(G), [~,G] = max(G); end
M = []; if isfield(IN,'M'), M = IN.M; end

if isfield(IN,'featind')
    featind = IN.featind;
else
    featind = true(1,size(Y,2));
end

if eIN || (~isfield(IN,'estimators') || isempty(IN.estimators) ) 

    if ~isfield(IN,'subgroup') || isempty(IN.subgroup)
        % Compute IN.beta from entire population
        [~,IN.estimators] = combat( Y(:,featind), G, M );
    else
        % Compute IN.beta from a subgroup of observations
        idxSubgroup = logical(IN.subgroup);
        [~,IN.estimators] = combat( Y(idxSubgroup,featind), G(idxSubgroup,:), M );
    end
end

Y = combat( Y(:,featind), G, M, IN.estimators);

Y=Y';

% =========================================================================
function [ Y, IN ] = FastICAObj( Y, IN )

if isempty(IN),eIN=true; else, eIN=false; end

if eIN|| ~isfield(IN,'G') || isempty(IN.G), error('No covariates defined in parameter structure'), end

if isfield(IN,'featind')
    featind = IN.featind;
else
    featind = true(1,size(Y,2));
end

if eIN || ~isfield(IN, 'sigICs') || isempty(IN.sigICs)

    if ~isfield(IN, 'subgroup') || isempty(IN.subgroup)
        [Y,IN] = cv_icaHarmonize(Y(IN.subgroup, featind), 1, IN);
    else
        [Y,IN] = cv_icaHarmonize(Y(:,featind), 1, IN);
    end
  
else

    Y = cv_icaHarmonize(Y(:,featind), 2, IN);
    
end