function [sY, IN, dT] = nk_PerfAdjForCovarsUsingPCAObj(Y, IN, S)
% =========================================================================
% function [adjT, IN] = nk_PerfAdjForCovarsUsingPCAObj(T, IN, S)
% =========================================================================
% ----- INPUTS -----
% Y :               Target matrix (row = patterns, cols = features), from 
%                   which the covariate effects will be removed 
% S :               Source Matrix (row = patterns, cols = features)
%                   Used to (1) determine which eigenvectors are correlated
%                   with IN.G, and (2) to reduce and reconstruct T
% IN.G :            Covariate matrix 
%
% Additional fields in IN include:
% IN.recon :        back-project adjusted target matrix to input space 
%                   (def: true)
% IN.varop :        operator for identification of correlated factors 
%                   (def: gt)
% IN.corrmeth :     identification method (1 = pearson (def), 
%                   2 = spearman, 3 = anova)
% IN.corrthresh :   cutoff for identification (def: 0.3)
% IN.DR :           Dimensionality reduction parameter structure (see
%                   nk_DimRed_config)
% IN.indX :         logical index vector to subcohort of IN.S
%
%
% ----- OUTPUTS -----
% adjT :            Adjusted target matrix 
% IN :              Parameter structure containing orig. and comp. params
% =========================================================================
% (c) Nikolaos Koutsouleris, 1/2016

% =========================== WRAPPER FUNCTION ============================ 
if ~exist('S','var'), S=[]; end
if iscell(Y) && exist('IN','var') && ~isempty(IN)
    sY = cell(1,numel(Y)); dT = [];
    for i=1:numel(Y), [sY{i}, IN] = PerfAdjForCovarsUsingPCAObj(Y{i}, IN, S); end
else
    [ sY, IN, dT ] = PerfAdjForCovarsUsingPCAObj( Y, IN, S );
end

end

% =========================================================================
function [adjT, IN, dT] = PerfAdjForCovarsUsingPCAObj(T, IN, S)
global VERBOSE

% Check existence of paramater structure
if ~exist('IN','var') || isempty(IN)          
    error('No parameter structure specified! Abort!');  
end
if ~isfield(IN,'recon') || isempty(IN.recon)
    IN.recon = true;
end
if ~isfield(IN,'varop') || isempty(IN.varop)
    IN.varop = 'gt';
end
    
compfl = false;
% Check whether variance removal parameters have been already computed
if ~isfield(IN,'ind0') || isempty(IN.ind0)
    
    compfl = true;
    % The source data matrix to compute the factorized matrix from
    if ~exist('S','var') || isempty(S)
        error('No training matrix specified in parameter structure'); 
    end
    % The source covariate matrix to compute correlation coefficients with
    % factorized data matrix
    if (~isfield(IN,'G') || isempty(IN.G))
        error('No target vector / matrix specified in parameter structure'); 
    end
    % The defaults correlation method
    if ~isfield(IN,'corrmeth') || isempty(IN.corrmeth)
        IN.corrmeth = 1;                                    
    end
    switch IN.corrmeth
        case 1
            corrmeth = 'pearson';
        case 2
            corrmeth = 'spearman';
        case 3
            corrmeth = 'anova';
    end

    % The default correlation strength 
    if ~isfield(IN,'corrthresh') || isempty(IN.corrthresh)
        IN.corrthresh = 0.3;                                
    end
    
    % if not otherwise specified use the entire source data matrix
    if ~isfield(IN,'indX') || isempty(IN.indX)
        IN.indX = true(size(S,1),1);
    end
    
    % if not otherwise specified use PCA in the 'percentage of
    % sum(eigenvalues) mode'

    % CV Note: this can potentially also be adjusted so that DR.RedMode is
    % directly set in _config
    if (~isfield(IN,'DR') || isempty(IN.DR) || ~isfield(IN.DR,'RedMode')) && (~isfield(IN, 'redmethod') || IN.redmethod == 1)
        IN.DR.DRsoft = 1; 
        IN.DR.RedMode = 'PCA';
        IN.DR.PercMode = IN.dimmode;  
    elseif (~isfield(IN,'DR') || isempty(IN.DR) || ~isfield(IN.DR,'RedMode')) && IN.redmethod == 2 
        IN.DR.RedMode = 'fastICA'; 
        IN.DR.PercMode = IN.dimmode; % unecessary, might affect nk_PerfRedObj
        IN.DR.DRsoft = 1; % unused, might affect nk_PerfRedPbj
    end
    
    % Run dimensionality reduction on the source matrix
    [dS, IN] = nk_PerfRedObj(S(IN.indX,:),IN);
    IN.C = zeros(size(dS,2),size(IN.G,2));
    IN.Pvalue = IN.C;

    % Identify correlated factors in the factorized source matrix
    switch corrmeth
        case {'pearson', 'spearman'}
            for i = 1:size(IN.G,2)
                % Determine correlations with covars in the source matrix
                [Cvalues,~,~,Pvalues] = nk_CorrMat(dS,IN.G(IN.indX,i),corrmeth);
                IN.C(:,i)=abs(Cvalues); IN.Pvalue(:,i)=Pvalues;
                IN.C(isinf(IN.C(:,i)) | isnan(IN.C(:,i)),i) = 0;
                if VERBOSE 
                    maxi = max(IN.C(:,i)); 
                    fprintf('\nWorking on covariate #%g => max=%g', i, maxi); 
                end
            end
        case 'anova'
            % This is the more powerful option if we have multiple covars
            % for which we would like to identify respective variance
            % components
            RES.X = [ones(size(IN.G(IN.indX,:),1),1) IN.G(IN.indX,:)];
            RES = nk_PerfANOVAObj(dS, RES);
            IN.C = sqrt(RES.R2);
            IN.Pvalue = RES.p;
            IN.C(isinf(IN.C) | isnan(IN.C)) = 0;
            if VERBOSE
                maxi = max(IN.C); 
                fprintf('\nDummy matrix => max=%g', maxi); 
            end
    end

    % Threshold eigenvariate correlations with covars
    switch IN.corrcrit
        case 'corr'
            IN.subthresh = any(single(feval(IN.varop, IN.C, IN.corrthresh)),2);
        case 'pval'
            IN.subthresh = any(single(feval(IN.varop, IN.Pvalue, IN.corrthresh)),2);
    end

    % Check whether correlated factors exist or not!
    if ~sum(IN.subthresh)
        IN.ind0 = true(1,size(dS,2));
        warning('\nNo variance components identified at %g threshold. Returning unchanged matrix!', IN.corrthresh);
        switch IN.recon
            case 1
              adjT = T;
            case 2
              if isempty(T) ||(isequal(T, S) && sum(IN.indX) == size(S,1))
                 adjT = dS; 
              else
                 adjT = nk_PerfRedObj(T,IN);
              end
        end
        dT = adjT;
        return
    end

    % Determine which eigenvariates have to be extracted
    IN.ind0 = ~IN.subthresh;
end
if ~any(IN.ind0) 
    warning('All eigenvariates meet threshold criterion [%s %g]!\nCheck your data and your settings.',IN.varop,IN.corrthresh);
end
    
% Now project target matrix to source matrix embedding space
if VERBOSE, fprintf(sprintf('\nProjecting target matrix to source %s space', IN.DR.RedMode)); end
if isempty(T) ||(exist('S','var') && isequal(T, S) && sum(IN.indX) == size(S,1))
    if exist('dS','var')
        dT = dS; 
    else
        % Run dimensionality reduction on the source matrix
        dT = nk_PerfRedObj(S(IN.indX,:),IN);
    end
else
    dT = nk_PerfRedObj(T,IN);
end

switch IN.recon
    case 1
        % Reconstruct target matrix only with / without extracted variance
        % componentsp
        if VERBOSE, fprintf('\nReconstructing target matrix based on selected components.'); end
        tT = bsxfun(@plus, IN.mpp.vec(:,IN.ind0)* dT(:,IN.ind0)' , IN.mpp.sampleMean')';
        adjT = zeros(size(dT,1),numel(IN.indNonRem)); adjT(:,IN.indNonRem) = tT;
    case 2
        if VERBOSE, fprintf('\nLimiting target matrix to selected components.'); end
        % Limited projected target matrix to selected variance components
        adjT = dT(:,IN.ind0);
end

if VERBOSE
    mx = IN.C(~IN.ind0,:); s0 = sum(~IN.ind0);
    if ~s0
        fprintf('???')
    end
    if compfl
        fprintf(['\nProcessing finished. %g/%g components with r %s %g (max: %g) were identified.', ...
            '\nRespective variance components was extracted from the input matrix.\n'],s0,size(IN.mpp.vec,2), IN.varop, IN.corrthresh, max(mx(:))); 
    else
        fprintf('\nProcessing finished. %g variance components were removed from input matrix', s0)
    end
end

end
