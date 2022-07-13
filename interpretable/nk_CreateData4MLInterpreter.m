function IN = nk_CreateData4MLInterpreter(RandFeats, MapIdx, Tr, Ts, covars, IN, nx)

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
        case {'median','medianflip','random'}
            IN.covars_rep = repmat(covars,IN.MLI.nperms,1);
    end
end

n = numel(MapIdx);

switch IN.MLI.method

    case 'posneg'
        % Determine extremes of the distribution
        upper = prctile(Tr(:,MapIdx), IN.MLI.upper_thresh);
        lower = prctile(Tr(:,MapIdx), IN.MLI.lower_thresh);
                
        if ~isinf(IN.MLI.nperms)
        
            P = repmat(Ts, IN.MLI.nperms, 1);
            N = repmat(Ts, IN.MLI.nperms, 1);
            I = false(IN.MLI.nperms, size(Tr,2));
            
            % Create modified instances of case
            for i=1:IN.MLI.nperms
                P(i, MapIdx(RandFeats(i,:))) = upper(RandFeats(i,:)); 
                N(i, MapIdx(RandFeats(i,:))) = lower(RandFeats(i,:));
                I(i, MapIdx(RandFeats(i,:))) = true;
            end
        else
            I = false(IN.MLI.max_iter, size(Tr,2)); iter = 1; completed = false;
            while iter <= IN.MLI.max_iter
                I(iter, MapIdx(RandFeats(iter,:))) = true;
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
        IN.X(nx).(Yocvstr){1} = P;
        IN.X(nx).(Yocvstr){2} = N;
        IN.X(nx).I = I;
      
    case {'median','medianflip','random'}
        switch IN.MLI.method
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
        if ~isinf(IN.MLI.nperms)
        
            M = repmat(Ts, IN.MLI.nperms, 1);
            I = false(IN.MLI.nperms, size(Tr,2));
            
            % Create modified instances of case
            for i=1:IN.MLI.nperms
                try
                    M(i, MapIdx(RandFeats(i,:))) = medi(RandFeats(i,:)); 
                    I(i, MapIdx(RandFeats(i,:))) = true;
                catch
                    fprintf('prob');
                end
            end
        else
            I = false(IN.MLI.max_iter, size(Tr,2)); iter = 1; completed = false;
            while iter <= IN.MLI.max_iter
                I(iter, MapIdx(RandFeats(iter,:))) = true;
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
        IN.X(nx).(Yocvstr) = M;
        IN.X(nx).I = I;
end

