function IN = nk_CreateData4MLInterpreter(Tr, Ts, covars, IN, nx)

n = size(Tr,2);

if IN.oocvflag
    Yocvstr = 'Yocv2';
else
    Yocvstr = 'Yocv';
end

if ~isempty(covars) && nx == 1
    switch IN.MLI.method
        case 'posneg'
            IN.covars_rep{1} = repmat(covars,IN.MLI.nperms,1);
            IN.covars_rep{2} = repmat(covars,IN.MLI.nperms,1);
        case 'median'
            IN.covars_rep = repmat(covars,IN.MLI.nperms,1);
    end
end
% If map is provided determine subspace for modification
if isfield(IN.MLI,'MAP') && IN.MLI.MAP.flag     
    cutoff = IN.MLI.MAP.cutoff;
    if isfield(IN.MLI.MAP,'percentile') && IN.MLI.MAP.percentmode
        cutoff = percentile(IN.MLI.MAP.map, IN.MLI.MAP.cutoff);
    end
    mapidx = return_imgind(IN.MLI.MAP.operator, cutoff, IN.MLI.MAP.map);
else
    mapidx = 1:n;
end

% Determine extremes of the distribution
n = numel(mapidx);
nfrac = ceil(n*IN.MLI.frac);

switch IN.MLI.method
    case 'posneg'

        upper = prctile(Tr(:,mapidx), IN.MLI.upper_thresh);
        lower = prctile(Tr(:,mapidx), IN.MLI.lower_thresh);
                
        if ~isinf(IN.MLI.nperms)
        
            P = repmat(Ts, IN.MLI.nperms, 1);
            N = repmat(Ts, IN.MLI.nperms, 1);
            I = false(IN.MLI.nperms, size(Tr,2));
            
            % Create modified instances of case
            for i=1:IN.MLI.nperms
                idx = randperm(n, nfrac);
                P(i, mapidx(idx)) = upper(idx); 
                N(i, mapidx(idx)) = lower(idx);
                I(i, mapidx(idx)) = true;
            end
        else
            I = false(IN.MLI.max_iter, size(Tr,2)); iter = 1; completed = false;
            while iter <= IN.MLI.max_iter
                I(iter, mapidx(randperm(n, nfrac))) = true;
                if sum( sum(I(1:iter,:)) >= IN.MLI.n_visited ) == n
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
            P = repmat(Ts, iter, 1);
            N = repmat(Ts, iter, 1);
            for i=1:iter
                P(i, I(i,:)) = upper(I(i,:)); 
                N(i, I(i,:)) = lower(I(i,:));
            end
        end
        % Remove duplicates
        [~, ix1] = unique(P-N,'rows'); 
        IN.X(nx).(Yocvstr){1} = P(ix1,:);
        IN.X(nx).(Yocvstr){2} = N(ix1,:);
        IN.X(nx).I = I(ix1,:);
        % Remove nans
        ix2 = isnan(IN.X(nx).(Yocvstr){1});
        IN.X(nx).(Yocvstr){1}(ix2) = 0;
        IN.X(nx).(Yocvstr){2}(ix2) = 0;

    case 'median'

        medi = prctile(Tr(:,mapidx), 50);
        if ~isinf(IN.MLI.nperms)
        
            M = repmat(Ts, IN.MLI.nperms, 1);
            I = false(IN.MLI.nperms, size(Tr,2));
            
            % Create modified instances of case
            for i=1:IN.MLI.nperms
                idx = randperm(n, nfrac);
                M(i, mapidx(idx)) = medi(idx); 
                I(i, mapidx(idx)) = true;
            end
        else
            I = false(IN.MLI.max_iter, size(Tr,2)); iter = 1; completed = false;
            while iter <= IN.MLI.max_iter
                I(iter, mapidx(randperm(n, nfrac))) = true;
                if sum( sum(I(1:iter,:)) >= IN.MLI.n_visited ) == n
                    completed = true;
                    break
                end
                iter=iter+1;
            end
            if ~completed
                warning('\n Not all features underwent the predefined amount of %g modifications after %g iterations.', IN.MLI.n_visited, IN.MLI.max_iter);
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
        [~, ix1] = unique(M,'rows');
        IN.X(nx).(Yocvstr){1} = M(ix1,:);
        IN.X(nx).I = I(ix1,:);
        % Remove nans
        ix2 = isnan(IN.X(nx).(Yocvstr){1});
        IN.X(nx).(Yocvstr){1}(ix2) = 0;
end

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