function IN = nk_CreateData4MLInterpreter(RandFeats, Tr, Ts, covars, IN, nx)

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
        case {'median','medianflip'}
            IN.covars_rep = repmat(covars,IN.MLI.nperms,1);
    end
end

n = numel(IN.MLI.MAP.mapidx);

switch IN.MLI.method
    case 'posneg'
        % Determine extremes of the distribution
        upper = prctile(Tr(:,IN.MLI.MAP.mapidx), IN.MLI.upper_thresh);
        lower = prctile(Tr(:,IN.MLI.MAP.mapidx), IN.MLI.lower_thresh);
                
        if ~isinf(IN.MLI.nperms)
        
            P = repmat(Ts, IN.MLI.nperms, 1);
            N = repmat(Ts, IN.MLI.nperms, 1);
            I = false(IN.MLI.nperms, size(Tr,2));
            
            % Create modified instances of case
            for i=1:IN.MLI.nperms
                P(i, IN.MLI.MAP.mapidx(RandFeats(i,:))) = upper(RandFeats(i,:)); 
                N(i, IN.MLI.MAP.mapidx(RandFeats(i,:))) = lower(RandFeats(i,:));
                I(i, IN.MLI.MAP.mapidx(RandFeats(i,:))) = true;
            end
        else
            I = false(IN.MLI.max_iter, size(Tr,2)); iter = 1; completed = false;
            while iter <= IN.MLI.max_iter
                I(iter, IN.MLI.MAP.mapidx(RandFeats(iter,:))) = true;
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
      
    case {'median','medianflip'}
        switch IN.MLI.method
            case 'median'
                medi = prctile(Tr(:,IN.MLI.MAP.mapidx), 50);
            case 'medianflip'
                centiles = nk_ComputePercentiles(Tr(:,IN.MLI.MAP.mapidx), Ts(:,IN.MLI.MAP.mapidx),'inverse');
                idxU = centiles>50;
                idxL = centiles<=50;
                centiles(idxU) = centiles(idxU) - 50;
                centiles(idxL) = centiles(idxL) + 50;
                medi = nk_ComputePercentiles(Tr(:,IN.MLI.MAP.mapidx), centiles, 'normal');
        end
        if ~isinf(IN.MLI.nperms)
        
            M = repmat(Ts, IN.MLI.nperms, 1);
            I = false(IN.MLI.nperms, size(Tr,2));
            
            % Create modified instances of case
            for i=1:IN.MLI.nperms
                M(i, IN.MLI.MAP.mapidx(RandFeats(i,:))) = medi(RandFeats(i,:)); 
                I(i, IN.MLI.MAP.mapidx(RandFeats(i,:))) = true;
            end
        else
            I = false(IN.MLI.max_iter, size(Tr,2)); iter = 1; completed = false;
            while iter <= IN.MLI.max_iter
                I(iter, IN.MLI.MAP.mapidx(RandFeats(iter,:))) = true;
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

