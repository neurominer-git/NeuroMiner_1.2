function IN = nk_CreateData4MLInterpreter(MLI, RandFeats, Tr, Ts, covars, IN, nx, curclass)

if IN.oocvflag
    Yocvstr = 'Yocv2';
else
    Yocvstr = 'Yocv';
end

method = MLI.method;
nperms = MLI.nperms;

if ~isempty(covars) && nx == 1
    switch method
        case 'posneg'
            IN.covars_rep{1} = repmat(covars,nperms,1);
            IN.covars_rep{2} = repmat(covars,nperms,1);
        case {'median','medianflip','random'}
            IN.covars_rep = repmat(covars,nperms,1);
    end
end

% MapIdx can be the full feature space or a selection based on the CVR and
% sign-based consistency signature chosen by the user. MapIdx is a logical
% vector indicating which variables should be used
MapIdx = MLI.Modality{nx}.MAP.mapidx{curclass}; 
fMapIdx = find(MapIdx);

n = sum(MapIdx);

switch method

    case 'posneg'
        % Determine extremes of the distribution
        upper = prctile(Tr(:,MapIdx), IN.MLI.upper_thresh);
        lower = prctile(Tr(:,MapIdx), IN.MLI.lower_thresh);
                
        if ~isinf(nperms)
        
            P = repmat(Ts, nperms, 1);
            N = repmat(Ts, nperms, 1);
            I = false(nperms, size(Tr,2));
            
            % Create modified instances of case
            for i=1:nperms
                
                % RandFeats is a random selection index to variables in MapIdx.
                % It can be selected by randomly indexing the original
                % input space variables or by using an atlas scheme that
                % groups input spaces variables together.
                % Select subsample of features by combining MapIdx and
                % current row of RandFeats. 

                rIdx = RandFeats(i,:);

                P(i, fMapIdx(rIdx)) = upper(rIdx); 
                N(i, fMapIdx(rIdx)) = lower(rIdx);
                I(i, fMapIdx(rIdx)) = true;
            end
        else
            I = false(MLI.max_iter, size(Tr,2)); iter = 1; completed = false;
            while iter <= MLI.max_iter
                I(iter, MapIdx(RandFeats(iter,:))) = true;
                if sum( sum(I(1:iter,:)) >= MLI.n_visited ) == n
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
        IN.X(nx).(Yocvstr){1} = P;
        IN.X(nx).(Yocvstr){2} = N;
        IN.X(nx).I = I;
      
    case {'median','medianflip','random'}
        switch method
            case 'median'
                medi = prctile(Tr(:,MapIdx), 50);
            case 'medianflip'
                centiles = nk_ComputePercentiles(Tr(:,MapIdx), Ts(:,MapIdx),'inverse');
                idxU = centiles>50;
                idxL = centiles<=50;
                centiles(idxU) = centiles(idxU) - 50;
                centiles(idxL) = centiles(idxL) + 50;
                medi = nk_ComputePercentiles(Tr(:,MapIdx), centiles, 'normal');
            case 'random'
                U = nk_CountUniques(Tr);
                medi = zeros(1,n);
                for i=1:n, medi(i) = U.UX{i}(randi(numel(U.UX{i}))); end
        end
        if ~isinf(nperms)
        
            M = repmat(Ts, nperms, 1);
            I = false(nperms, size(Tr,2));
            
            % Create modified instances of case
            for i=1:nperms
                rIdx = RandFeats(i,:);
                M(i, fMapIdx(rIdx)) = medi(rIdx); 
                I(i, fMapIdx(rIdx)) = true;
            end
        else
            I = false( MLI.max_iter, size(Tr,2) ); iter = 1; completed = false;
            while iter <= MLI.max_iter
                I(iter, MapIdx(RandFeats(iter,:))) = true;
                if sum( sum(I(1:iter,:)) >= MLI.n_visited ) == n
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
        IN.X(nx).(Yocvstr) = M;
        IN.X(nx).I = I;
end