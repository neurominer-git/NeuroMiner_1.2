function [M, I] = nkdev_CreateData4ModelInterpreterMedian(Tr, Ts, IN, Px)

n = size(Tr,2);

% If map is provided determine subspace for modification
if isfield(IN,'MAP') && ~isempty(IN.MAP)
    
    IN.MAP.map = return_inputspace(IN.MAP.map, IN, Px);
    cutoff = IN.MAP.cutoff;
    if isfield(IN.MAP,'percentmode') && IN.MAP.percentmode
        cutoff = percentile(IN.MAP.map, IN.MAP.cutoff);
    end
    mapidx = return_imgind(IN.MAP.operator, cutoff, IN.MAP.map);
else
    mapidx = 1:n;
end

% Determine extremes of the distribution
n = numel(mapidx);
medi = prctile(Tr(:,mapidx), 50);
%medi = zeros(1,numel(mapidx));

nfrac = ceil(n*IN.frac);

if ~isinf(IN.nperms)

    M = repmat(Ts, IN.nperms, 1);
    I = false(IN.nperms, size(Tr,2));
    
    % Create modified instances of case
    for i=1:IN.nperms
        idx = randperm(n, nfrac);
        M(i, mapidx(idx)) = medi(idx); 
        I(i, mapidx(idx)) = true;
    end
else
    I = false(IN.max_iter, size(Tr,2)); iter = 1; completed = false;
    while iter <= IN.max_iter
        I(iter, mapidx(randperm(n, nfrac))) = true;
        if sum( sum(I(1:iter,:)) >= IN.n_visited ) == n
            completed = true;
            break
        end
        iter=iter+1;
    end
    if ~completed
        warning('\n Not all features underwent the predefined amount of %g modifications after %g iterations.', IN.n_visited, IN.max_iter);
    else
        I(iter,:)=[];
        iter=iter-1;
    end
    M = repmat(Ts, iter, 1);
    for i=1:iter
        M(i, I(i,:)) = medi(I(i,:)); 
    end
end
% Remove duplicates
[~, ix] = unique(M,'rows');
M = M(ix,:);
I = I(ix,:);
ix = isnan(M);
M(ix) = 0;

function imgind = return_imgind(typthresh, thresh, img)

if length(thresh) > 1
    switch typthresh
        case 1
            imgind = (img < thresh(1) | img > thresh(2)); 
        case 2
            imgind = (img <= thresh(1) | img >= thresh(2)); 
        case 3
            imgind = (img > thresh(1) | img < thresh(2)); 
        case 4
            imgind = (img >= thresh(1) | img <= thresh(2)); 
    end
else
    switch typthresh
        case 1
            imgind = img < thresh; 
        case 2
            imgind = img <= thresh;
        case 3
            imgind = img > thresh;
        case 4
            imgind = img >= thresh;
        case 5
            imgind = img == thresh;
    end
end
imgind = find(imgind);
% -------------------------------------------------------------------------
function nmW = return_inputspace(nmW, IN, naPX)

if isfield(IN,'DRMODE') && ~isempty(IN.DRMODE)
    
    switch IN.DRMODE.DR.RedMode                                        
        case {'PCA','RobPCA','SparsePCA'}
            if strcmp(IN.DRMODE.DR.RedMode,'RobPCA')
                DRsoft = 1; 
            else
                DRsoft = IN.DRMODE.DR.DRsoft;
            end
            % Project back to input space
            switch DRsoft
                case 0
                    nmW = reconstruct_data(nmW, naPX.mpp);
                case 1
                    nmW = nmW * naPX.mpp.vec';     
            end
        case {'optNMF','NeNMF','NNMF','PLS','LPP', 'NPE', 'LLTSA', ...
                'SPCA', 'PPCA', 'FA', 'FactorAnalysis', 'NCA', 'MCML', 'LMNN'}
             nmW = nmW * naPX.mpp.vec';
        otherwise
             error('Reconstruction of data is not supported for this technique.');
    end
    
    % If features had been removed prior to dimensionality
    % reduction take care that you recover the original space
    if isfield(naPX,'indNonRem') && ~isempty(naPX.indNonRem) && sum(~naPX.indNonRem) > 0
        tmW = zeros(size(naPX.indNonRem')); 
        tmW(naPX.indNonRem) = nmW; nmW = tmW; 
    end

end