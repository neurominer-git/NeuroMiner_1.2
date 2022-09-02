function [sY, IN] = nk_PerfImputeObj(Y, IN)
% =========================================================================
% FORMAT function [sY, IN] = nk_PerfImputeObj(Y, IN)
% =========================================================================
% Performes imputation using either single-subject median replacement, 
% feature-wise mean replacement or multivariate distance-based NN median 
% imputation. If you want to use hamming or euclidean distances you should 
% scale, unit-normalize or standardize the data first, otherwise the
% distance measure will be dominated by high-variance features.
% 
% Input/Output arguments:
% IN.blockind       : [ ]            Feature block index vector (boolean).
%                                    if [] all features in matrix are used
% IN.method         : [ singlemean ] single-subject median replacement if
%                                    NaN value in feature block
%                     [ mean ]       NaN value is replaced by mean of value
%                                    within the given feature
%                     [ euclidean, cityblock, seuclidean, cosine, ... ] 
%                                   NN-based replacement of NaN using
%                                   the median of the IN.K most similar in-
%                                   stances with finite values for the given 
%                                   NaN value. Similarity is determined 
%                                   from IN.X using IN.method. Requires 
%                                   pdist2 which is available through 
%                                   the MATLAB statistics toolbox  
% sY                :               The imputed data matrix
% Changes:
% 25/03/2022        : Improved the selection of training subjects for
%                     imputation by sorting the training matrix and
%                     identifying observation without NaNs. 
% 28/07/2022        : Further changes introduced to avoid imputation errors
%                     when no columns can be found without missings to 
%                     establish similarity. Then, cases are iteratively 
%                     removed from the imputation learning sample until the 
%                     criterion minnumcols (currently fixed at 0.5) is 
%                     fulfilled.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 08/2022

% =========================== WRAPPER FUNCTION ============================ 
if iscell(Y) 
    if ~exist('IN','var'), IN.X = Y{1}; end
    sY = cell(1,numel(Y)); 
    for i=1:numel(Y), [sY{i}, IN] =  PerfImputeObj(Y{i}, IN); end
else
    if ~exist('IN','var'), IN.X = Y; end
    [ sY, IN ] = PerfImputeObj(Y, IN );
end
% =========================================================================
function [sY, IN] = PerfImputeObj(Y, IN)
global VERBOSE 

% Check for cases that are completely NaN and remove them temporarily
[Y, ~, IxNaN] = nk_ManageNanCases(Y, []);

%% Prepare
[m, n] = size(Y);

if ~isfield(IN,'method') || isempty(IN.method), IN.method = 'euclidean'; end
if ~isfield(IN,'k') || isempty(IN.k), IN.k = 7; end
if ~isfield(IN,'blockind') || isempty(IN.blockind), IN.blockind = 1:n; end
if ~isfield(IN,'minnumcols') || isempty(IN.minnumcols), IN.minnumcols = 0.75; end
if ~isfield(IN,'X') && ~strcmp(IN.method,'singlemean'), ...
        error('The training data matrix is missing from the input parameters!'); end
sY      = Y;
tY      = Y(:,IN.blockind);
stY     = tY;
indnan  = isnan(tY);
snan    = sum(indnan,2)==0;
minnumcols = ceil(n * IN.minnumcols);
if VERBOSE, fprintf('\tImpute missing values');end
ll=0; 
switch IN.method
    case 'singlemean'
        mn = nm_nanmedian(tY,2); 
        for i = 1:m
            if snan(i), if VERBOSE,fprintf('.'); end; continue; end
            if VERBOSE,fprintf('+'),  end
            stY(i,indnan(i,:)) = mn(i);
            ll=ll+1;
        end
    case 'mean'
        tX = IN.X(:,IN.blockind); 
        mn = nm_nanmean(tX); 
        for i = 1:m
            if snan(i), if VERBOSE,fprintf('.'); end; continue; end
            if VERBOSE,fprintf('+'),  end
            stY(i,indnan(i,:)) = mn(indnan(i,:));
            ll=ll+1;
        end
    case {'euclidean','cityblock','seuclidean','cosine','mahalanobis','jaccard','hamming','hybrid'}
        tX = IN.X(:,IN.blockind);
        IN.C = nan(m,1);
        if strcmp(IN.method,'hybrid')
            R = nk_CountUniques(tX);
            indNom = (R.U <= IN.hybrid.cutoff)';
        end
        for i = 1:m
            if snan(i), continue; end
            if VERBOSE; fprintf('.'); end
            % Find which features are NaN in the given case i
            ind_Yi = find(indnan(i,:)); 
            n_ind_Yi = numel(ind_Yi);

            for j=1:n_ind_Yi
                % Get training cases which do not have NaNs in the given
                % column 
                indnan_Xi = find(~isnan(tX(:,ind_Yi(j))));
                if ~sum(indnan_Xi)
                    fprintf('\n');warning('I did not find observations with non-missing values for feature %g. Check your settings!', j)
                end
                indi_Yi = 1:size(tX,2); indi_Yi(ind_Yi(j))=[];
                
                % Find cases without missings 
                tXX = tX(indnan_Xi, indi_Yi);
                indnan_Xj = isnan(tXX);
                [~, ind_c] = sort(sum(indnan_Xj),'ascend');
                [~, ind_r] = sort(sum(isnan(tXX), 2),'ascend');
                idx_c = sum(isnan(tX(indnan_Xi, indi_Yi(ind_c))))==0;
                while sum(idx_c) < minnumcols
                    % if there is no column without missings start removing
                    % cases from the imputation training sample
                    ind_r(end)=[];
                    idx_c = sum(isnan(tX(indnan_Xi(ind_r), indi_Yi(ind_c))))==0;
                    if numel(indnan_Xi)<IN.k
                        warning('Remaining samples in training data < k ');
                    end
                end
                indi_Yi = indi_Yi(ind_c(idx_c));
                indnan_Xi = indnan_Xi(ind_r);
                Xj = tX(indnan_Xi, indi_Yi);
    
                % Compute distance metric
                switch IN.method
                    case 'seuclidean'
                        S = nm_nanstd(Xj); S(S==0) = min(S(S~=0));
                        D = pdist2(Xj, tY(i,indi_Yi), 'seuclidean',S)';
                    case 'mahalanobis'
                        %C = nancov(Xj(:,indi_Yi)); C(C==0) = min(C(C~=0));
                        D = pdist2(Xj, tY(i,indi_Yi), 'mahalanobis')';
                    case 'hybrid'
                        % Identify nominal features using predefined
                        % variable
                        % Compute distances in nominal features
                        indZ1 = indi_Yi & indNom;
                        indZ2 = indi_Yi & ~indNom; 
                        D1 =  pdist2(tX(indnan_Xi, indZ1), tY(i,indZ1),IN.hybrid.method1);
                        % Compute distances in ordinal / continuous features;
                        D2 =  pdist2(tX(indnan_Xi, indZ2), tY(i,indZ2),IN.hybrid.method2);
                        D = mean([nk_PerfScaleObj(D1) nk_PerfScaleObj(D2)],2);
                    otherwise
                        D = pdist2(Xj, tY(i,indi_Yi),IN.method)';
                end
                % Sort training cases according to their proximity to the
                % test case whose value will be imputed.
                [Ds, ind] = sort(D,'ascend');
                if numel(Ds) < IN.k, kx = numel(Ds); else, kx = IN.k; end
                % Here we take the unweighted median of the kx imputation
                % training cases
                mn = median( tX(indnan_Xi(ind(1:kx)), ind_Yi(j)) );
                if isnan(mn)
                    fprintf('\n');warning('NaNs remain in the the current observation. Check your settings!')
                end
                stY(i,ind_Yi(j)) = mn;
                %if VERBOSE,fprintf('+'), end
                ll=ll+1; 
            end
            if isfield(IN,'compute_corr') && IN.compute_corr
                IN.C(i) = mean(nk_CorrMat(stY(i,indi_Yi)',Xj(ind(1:kx),indi_Yi)'));
            end
        end
end
if VERBOSE, fprintf('\n\t\t\t%g subject(s) with NaNs, total of %g NaN replaced.',sum(~snan),ll); end
sY(:,IN.blockind) = stY;
if sum(isnan(sY(:))) 
    warning(sprintf('\nNot all NaNs imputed!'))
end
% If you remove completely NaN cases temporarily add them back now to the imputed data.
sY = nk_ManageNanCases(sY, [], IxNaN);