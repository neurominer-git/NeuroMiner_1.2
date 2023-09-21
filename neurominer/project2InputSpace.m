function  bpX = project2InputSpace(X, PREPROC, prevP)
nPREPROC = PREPROC;
nPX = prevP;
nA = length(nPX);

%%%%%%%%%%% PROJECT DATA BACK TO INPUT SPACE %%%%%%%%%%

% Revert dimensionality reduction if previously used
% Find Dimensionality reduction parameters
reducedimfl = false;

for a = nA:-1:1
   
    % Adjust pnt according to parameter combination in current
    % preprocessing step

    if iscell(nPX{a})
        naPX = nPX{a}{1};
    else
        naPX = nPX{a};
    end

    switch nPREPROC.ACTPARAM{a}.cmd

        case 'scale'
            % Very frequently users opt to scale features after
            % factorization/dimensionality reduction. In these cases
            % scaling has to be reverted prior to back-projection
            if ~reducedimfl %&& decompfl(1)
                IN = nPREPROC.ACTPARAM{a}.SCALE;
                IN.minY = naPX.minY; IN.maxY = naPX.maxY;
                IN.revertflag = true;
                nmW(nmW==0)=NaN; nmW = nk_PerfScaleObj(nmW', IN)'; nmW(isnan(nmW))=0;
            end

        case {'reducedim','remvarcomp'}

            if isfield(naPX,'recon') && naPX.recon==1
                fprintf('-');
            else
                if isfield(naPX.mpp,'vec')
                    redvec = naPX.mpp.vec;
                elseif isfield(naPX.mpp,'factors')
                    redvec = naPX.mpp.factors{1};
                elseif isfield(naPX.mpp,'u')
                    redvec = naPX.mpp.u;
                elseif isfield(naPX.mpp,'M')
                    redvec = naPX.mpp.M;
                elseif isfield(naPX.mpp,'W')
                    redvec = naPX.mpp.W;
                elseif isfield(naPX.mpp,'network')
                    error('Autoencoder reconstructions not supported!')
                end

                if isfield(naPX,'ind0')
                    ind0 = naPX.ind0;
                    DR = naPX.DR;
                else
                    ind0 = 1:size(redvec,2);
                    DR = nPREPROC.ACTPARAM{a}.DR;
                end

                mpp.vec = redvec(:,ind0);

                switch DR.RedMode
                    case {'PCA','RobPCA','SparsePCA'}
                        if strcmp(DR.RedMode,'RobPCA')
                            DRsoft = 1;
                        else
                            DRsoft = DR.DRsoft;
                        end
                        % Project back to input space
                        switch DRsoft
                            case 0
                                nmW = reconstruct_data(nmW, mpp);
                                nmP = logical(reconstruct_data(nmP, mpp));
                                if ~isempty(PA), nmPA = reconstruct_data(nmPA, mpp); end
                            case 1
                                nmW = mpp.vec * nmW;
                                nmP = logical(mpp.vec * nmP);
                                if ~isempty(PA), nmPA = mpp.vec * nmPA; end
                        end
                    case {'optNMF','NeNMF','NNMF','PLS','LPP', 'NPE', 'LLTSA', 'SPCA', 'PPCA', 'FA', 'FactorAnalysis', 'NCA', 'MCML', 'LMNN'}
                        nmW = mpp.vec * nmW;
                        nmP = logical(mpp.vec * nmP);
                        if ~isempty(PA), nmPA = mpp.vec * nmPA; end

                    otherwise
                        error('Reconstruction of data is not supported for this technique.');
                end

                % If features had been removed prior to dimensionality
                % reduction take care that you recover the original space
                if isfield(naPX,'indNonRem') && ~isempty(naPX.indNonRem) && sum(~naPX.indNonRem) > 0
                    tmW = zeros(size(naPX.indNonRem')); tmP = tmW;
                    tmW(naPX.indNonRem) = nmW; nmW = tmW; tmP(naPX.indNonRem) = nmP; nmP = tmP;
                    if ~isempty(PA), tmPA = tmW; tmPA(naPX.indNonRem) = nmPA; nmPA = tmPA; end
                    % Don't forget to adjust the feature masks and the
                    % indices to modalities in case of fused feature spaces
                    %tlFuVI = false(length(naPX.indNonRem),1); tlFuVI(naPX.indNonRem) = lFuVI; lFuVI = tlFuVI;
                    %tlVI = true(length(naPX.indNonRem),1); tlVI(naPX.indNonRem) = fVI; fVI = tlVI;
                    clear tmW tmP %tlFuVI tlVI;
                end
                reducedimfl = true;
            end

        case {'elimzero','extfeat','extdim'}
            % Find eliminated feature vector
            if isfield(naPX,'NonPruneVec')
                IND = 'NonPruneVec';
            elseif isfield(naPX,'indNonRem')
                IND = 'indNonRem';
            else
                IND = 'ind';
            end
            % Retrieve index vector to eliminated features
            if size(naPX.(IND),2)>1 && size(naPX.(IND),1)>1
                pIND = naPX.(IND)(:,inp.curlabel);
            else
                pIND = naPX.(IND);
            end
            % Eliminated features have to be re-introduced
            tmW = zeros(length(pIND),1); tmW(pIND) = nmW; nmW = tmW;
            tmP = false(length(pIND),1); tmP(pIND) = nmP; nmP = tmP;
            if ~isempty(PA)
                tmPA = tmW; tmPA(naPX.(IND)) = nmPA; nmPA = tmPA;
            end
            
    end
end
bpX = nmPA;
end