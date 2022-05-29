function [Results, FileNames, RootPath] = nk_MLInterpreter(inp)
% =========================================================================
% FORMAT Results = nk_MLInterpreter(inp)
% =========================================================================
% Interpretable ML module
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, last modified 08/2020

global SVM RFE MODEFL CV SCALE SAV CVPOS MLI 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
FullPartFlag    = RFE.ClassRetrain;
switch inp.analmode
    case 0
        ovrwrt  = inp.ovrwrt;                       % overwrite existing data files
    case 1
        matfiles = inp.matfiles;                    % MLI files
end
multiflag       = false; if strcmp(MODEFL,'classification'), multiflag = inp.multiflag; end
saveparam       = inp.saveparam;
loadparam       = inp.loadparam;
nclass          = inp.nclass;
analysis        = inp.analysis;
GridAct         = inp.GridAct;
batchflag       = inp.batchflag;
algostr         = GetMLType(SVM);

% Setup CV2 and CV1 counters and data containers:
[ix, jx]        = size(CV.TrainInd);
[iy, jy]        = size(CV.cvin{1,1}.TrainInd);
%interpAll       = size(inp)
ll              = 1; 
totLearn        = 0;
oocvflag        = inp.oocvflag;
inp.MLI         = MLI;

if ~exist('GridAct','var') || isempty(GridAct), GridAct = nk_CVGridSelector(ix,jx); end
if ~exist('batchflag','var') || isempty(batchflag), batchflag = false; end

% Check and transform labels if needed
inp = nk_ApplyLabelTransform( SCALE, MODEFL, inp );

% Check whether you have to do prediction detrending for regression models
detrendfl = false;
if isfield(SVM,'Post') && isfield(SVM.Post,'Detrend') && SVM.Post.Detrend && strcmp(MODEFL,'regression')
    detrendfl = true;
end

% Check whether you have to perform label imputation and set flags
IMPUTE.flag = false;
if iscell(inp.PREPROC), iPREPROC = inp.PREPROC{1}; else, iPREPROC = inp.PREPROC; end    
if isfield(iPREPROC,'LABELMOD') && isfield(iPREPROC.LABELMOD,'LABELIMPUTE')
    IMPUTE = iPREPROC.LABELMOD.LABELIMPUTE; 
    IMPUTE.flag = true; 
end

BINMOD = iPREPROC.BINMOD;

CVPOS.fFull = FullPartFlag;

FileNames = cell(ix,jx); fnd=false;
nM = numel(inp.X);
RandFeats = struct('I',[]);
for h=1:nclass
    switch MODEFL
        case 'classification'
            inp.MLI.RangePred(h) = range(analysis.BinClass{h}.mean_predictions);
        case 'regression'
            inp.MLI.RangePred = range(analysis.Regr.mean_predictions);
    end
    
    for nx=1:nM
        nY = size(inp.X(nx).Y,2);
        if inp.oocvflag
            [mY2,nY2] = size(inp.X(nx).Yocv);
        else
            [mY2,nY2] = size(inp.X(nx).Y);
        end
        switch MODEFL
            case 'classification'
                Results.BinResults(h).Modality(nx).Y_mapped = zeros(mY2,nY2,ix);
                Results.BinResults(h).Modality(nx).Y_mapped_ciu = zeros(mY2,nY2,ix);
                Results.BinResults(h).Modality(nx).Y_mapped_cil = zeros(mY2,nY2,ix);
                Results.BinResults(h).Modality(nx).Y_mapped_std = zeros(mY2,nY2,ix);
            case 'regression'
                Results.RegrResults.Modality(nx).Y_mapped = zeros(mY2,nY2,ix);
                Results.RegrResults.Modality(nx).Y_mapped_ciu = zeros(mY2,nY2,ix);
                Results.RegrResults.Modality(nx).Y_mapped_cil = zeros(mY2,nY2,ix);
                Results.RegrResults.Modality(nx).Y_mapped_std = zeros(mY2,nY2,ix);
        end
        
        % If map is provided determine subspace for modification
        if isfield(inp.MLI,'MAP') && inp.MLI.MAP.flag && isfield(inp,'visdata')
            maptype = inp.MLI.MAP.map; inp.MLI.MAP.map=[];
            switch maptype
                case 'cvr'
                    inp.MLI.MAP.map{h} = inp.visdata{inp.curmodal,inp.curlabel}.CVRatio{h};
                case 'p_sgn'
                    inp.MLI.MAP.map{h} = inp.visdata{inp.curmodal,inp.curlabel}.SignBased_CV2_p_uncorr{h};
                case 'p_FDR_sgn'
                    inp.MLI.MAP.map{h} = inp.visdata{inp.curmodal,inp.curlabel}.SignBased_CV2_p_fdr{h};
            end
            cutoff = inp.MLI.MAP.cutoff;
            if isfield(inp.MLI.MAP,'percentile') && inp.MLI.MAP.percentmode
                cutoff = prctile(inp.MLI.MAP.map{h}, inp.MLI.MAP.cutoff);
            end
            inp.MLI.MAP.mapidx = return_imgind(inp.MLI.MAP.operator, cutoff, inp.MLI.MAP.map{h});
        else
            inp.MLI.MAP.mapidx = 1:nY;
        end
    end
end

nYmap = numel(inp.MLI.MAP.mapidx);
nfrac = ceil(nYmap*inp.MLI.frac);
permfile = fullfile(inp.rootdir,[SAV.matname '_MLIpermmat_ID' inp.id '.mat']);
if exist(permfile,'file') && inp.ovrwrtperm == 2
    fprintf('\nLoading %s', permfile);
    load(permfile)
else
    clc
    sprintf('Generating random feature subspaces ... \n\n\n');
    for h=1:nclass
        if strcmp(MODEFL,'classification') && nclass >1, fprintf('\n\tBinary classifier #%g', h); end
        for nx=1:nM
            if nM > 1, fprintf('\n\t\tModality #%g', nx); end
            [ RandFeats(h, nx).I, inp.MLI.nperms ] = uperms( inp.MLI.MAP.mapidx, inp.MLI.nperms, nfrac);   
        end
    end 
    fprintf('\nSaving to %s', permfile);
    save(permfile, 'RandFeats');
end
ol=1;
% =========================================================================
for f=1:ix % Loop through CV2 permutations

    for d=1:jx % Loop through CV2 folds
        
        fprintf('\n--------------------------------------------------------------------------')
        if ~GridAct(f,d) 
            ll=ll+1;
            fprintf('\nSkipping CV2 [%g,%g] (user-defined).',f,d)
            continue 
        end
        if ~fnd
            ll_start=ll; fnd=true; 
        end
        % Prepare variables for current CV2 partition
        CVPOS.CV2p = f;
        CVPOS.CV2f = d;
        if oocvflag
            tInd =  1:numel(inp.nOOCVsubj);
        else
            tInd = CV.TestInd{f,d};
        end
        inp.f = f; inp.d = d; 
        operm = f; ofold = d;
        inp.ll = ll;

        % Create OOCV partition file path
        oMLIpath = nk_GenerateNMFilePath(inp.rootdir, SAV.matname, inp.datatype, inp.multlabelstr, inp.varstr, inp.id, operm, ofold);
        OptModelPath = nk_GenerateNMFilePath( inp.saveoptdir, SAV.matname, 'OptModel', [], inp.varstr, inp.id, operm, ofold);

        switch inp.analmode
            case 0
                 % Compute from scratch
                 % Now, check whether file exist and can be loaded
                 loadfl = false;
                 if exist(oMLIpath,'file') && ~ovrwrt && ~batchflag
                    
                    [~, onam] = fileparts(oMLIpath);
                    fprintf('\nMLIdatamat found for CV2 [%g,%g]:',f,d)
                    fprintf('\nLoading: %s',onam)
                    try
                        load(oMLIpath)
                        loadfl = true;
                    catch
                        fprintf('\nCould not open file. May be corrupt. Recompute CV2 partition [%g,%g].',f,d);
                        loadfl = false;
                    end
                    
                elseif exist(oMLIpath,'file') && batchflag
                     % in batch mode we do not compute statistics across the
                    % CV2 partitions
                    [~, onam] = fileparts(oMLIpath);
                    fprintf('\nINTERdatamat found for CV2 [%g,%g]:\n%s',f,d,onam)
                    fprintf('\nBatch mode detected. Continue.')
                    [RootPath, FileNames{f,d}] = fileparts(oMLIpath); 
                    load(oMLIpath)
                    loadfl=true;
                end
                   
                if ~loadfl 

                    % Parameter flag structure for preprocessing
                    paramfl = struct('use_exist',inp.loadparam, ...
                                     'found', false, ...
                                     'write', inp.saveparam, ...
                                     'writeCV1', inp.saveCV1, ...
                                     'multiflag', multiflag);
                    
                    % Perform preprocessing of CV1/CV2 data
                    [ inp, contfl, analysis, mapY, GD, MD, Param, paramfl, mapYocv ] = nk_ApplyTrainedPreproc(analysis, inp, paramfl);

                    if contfl, continue; end

                    % Can we use a pretrained model saved on disk?
                    fndMD = false; 
                    if loadparam && isfield(inp,'optmodelmat') && exist(inp.optmodelmat{operm,ofold},'file')
                        fprintf('\nLoading OptModel: %s', inp.optmodelmat{operm,ofold});
                        load(inp.optmodelmat{operm,ofold},'MD'); fndMD = true; 
                    end
                    if ~fndMD, MD = cell(nclass,1); end
                    predOrig = cell(nclass,1);
                    % -----------------------------------------------------------------
                    for h=1:nclass  % Loop through binary comparisons
    
                        if nclass > 1, fprintf('\n\n*** %s #%g ***',algostr, h); end

                        %% Step 1: Get optimal model parameters
                        % Retrieve optimal parameters from precomputed analysis structure
                        % Differentiate according to binary or multi-group mode
                        [Ps, Pspos, nP, Pdesc] = nk_GetModelParams2(analysis, multiflag, ll, h, inp.curlabel);

                        if ~fndMD , MD{h} = cell(nP,1); end
                        
                        for m = 1 : nP % Loop through parameter combinations
                            
                            % Create model parameter string
                            cPs = Ps(m,:); sPs = nk_PrepMLParams(Ps, Pdesc, m);                            
                            P_str = nk_DefineMLParamStr(cPs, analysis.Model.ParamDesc, h);
                            
                            %% Step 2: Apply trained model to original target data 
                            % Optionally, retrain every base learner in current CV1
                            % [k,l] partition (across different grid positions, if available)
                            if ~fndMD,MD{h}{m} = cell(iy,jy); end

                            for k=1:iy % Loop through CV1 permutations

                                for l=1:jy % Loop through CV1 folds
                                    
                                    CVPOS.CV1p = k;
                                    CVPOS.CV1f = l;

                                    % Get feature feature subspace mask for current 
                                    % parameter grid position
                                    Fkl = GD.FEAT{Pspos(m)}{k,l,h};

                                    % Determine number of features in mask and
                                    % convert feature mask to logical index, if needed
                                    ul=size(Fkl,2); totLearn = totLearn + ul;
                                    if ~islogical(Fkl), F = Fkl ~= 0; else, F = Fkl; end

                                    % Get data pointers for current dichotomization
                                    TrInd = mapY.TrInd{k,l}{h};
                                    CVInd = mapY.CVInd{k,l}{h};
                                   
                                    % Set the pointer to the correct mapY 
                                    % (preprocessed data) shelf
                                    for n=1:numel(paramfl)
                                        pnt = 1;
                                        if ~BINMOD
                                             if isfield(paramfl{n},'PREPROC') && ...
                                               isfield(paramfl{n},'PXfull') && ...
                                               ~isempty(paramfl{n}.P{1})
                                                pnt = m;
                                                break   
                                            end
                                        else
                                            if isfield(paramfl{n},'PREPROC') && ...
                                               isfield(paramfl{n},'PXfull') && ...
                                               ~isempty(paramfl{n}.P{h})
                                                %Actual parameter node
                                                pnt = m; 
                                                break   
                                            end
                                        end
                                    end

                                    % get training data using pointers
                                    % Either (a) only CV1 training data, (b) CV1
                                    % training and test data, (c) CV1 training & test
                                    % data as well as CV2 test data              
                                    if BINMOD, hix = h; else, hix = 1; end
                                    if inp.oocvflag
                                        [ TR , CV1, CV2, OCV ] = nk_ReturnAtOptPos(mapY.Tr{k,l}{hix},  mapY.CV{k,l}{hix}, mapY.Ts{k,l}{hix}, mapYocv.Ts{k,l}{hix}, Param{1}(k,l,hix), pnt);
                                        uD = zeros(size(OCV,1),ul);
                                    else 
                                        [ TR , CV1, CV2 ] = nk_ReturnAtOptPos(mapY.Tr{k,l}{hix},  mapY.CV{k,l}{hix}, mapY.Ts{k,l}{hix}, [], Param{1}(k,l,hix), pnt);
                                        uD = zeros(size(CV2,1),ul);
                                    end
                                    if FullPartFlag, TR = [ TR; CV1]; end
 
                                   % Get and build label info
                                    modelTrL = mapY.TrL{k,l}{h};                         
                                    if FullPartFlag 
                                        modelTrL = [modelTrL; mapY.CVL{k,l}{h}]; 
                                        TrInd = [TrInd; CVInd]; 
                                    end

                                    % Loop through feature subspaces
                                    if ~fndMD 
                                        MD{h}{m}{k,l} = cell(ul,1); 
                                        fprintf(['\nRetrain models in CV2 [%g,%g], ' ...
                                            'CV1 [%g,%g], %g %s, (total # learners: %g), ML params [%s] ==> '], ...
                                            f, d, k, l, ul, algostr, totLearn, P_str)
                                        % Impute labels if needed
                                        [modelTrL, TR, TrInd] = nk_LabelImputer( modelTrL, TR, TrInd, sPs, IMPUTE);
                                        TR = TR(TrInd,:);
                                    else
                                        fprintf(['\nUse precomputed models in CV2 [%g,%g], ' ...
                                            'CV1 [%g,%g], %g %s, (total # learners: %g), ML params [%s] ==> '], ...
                                            f, d, k, l, ul, algostr, totLearn, P_str)
                                    end

                                    % Loop through feature subspaces
                                    for u=1:ul
                                        % Extract features according to mask
                                        TR_star = nk_ExtractFeatures(TR, F, [], u);
                                        
                                        if ~fndMD
                                            fprintf('Compute optimal model...');
                                            [~, MD{h}{m}{k,l}{u}] = nk_GetParam2(TR_star, modelTrL, sPs, 1);
                                            fprintf(' done.');
                                        end

                                        fprintf(' Predict CV2 test data...')
                                        % Apply model to CV2 test data 
                                        if inp.oocvflag
                                            [~, ~, uD(:,u)] = nk_GetTestPerf(OCV, ones(size(OCV,1),1), F(:,u), MD{h}{m}{k,l}{u}, TR_star, 1);
                                        else
                                            [~, ~, uD(:,u)] = nk_GetTestPerf(CV2, ones(size(CV2,1),1), F(:,u), MD{h}{m}{k,l}{u}, TR_star, 1);

                                        end
                                        fprintf(' done.')
                                    end
                                    
                                    if detrendfl 
                                        fprintf(' Detrend CV2 predictions.')
                                        beta = GD.Detrend{Pspos(m)}.beta;
                                        p = GD.Detrend{Pspos(m)}.p;
                                        uD(:,u) = nk_DetrendPredictions2(beta, p, uD(:,u)); 
                                    end    
                                    predOrig{h} = [predOrig{h} uD];
                                end
                            end
                        end
                    end
                    
                    %% Step 3: Create artificial versions of target cases
                    % Prepare arrays
                    nTs = numel(tInd);
                  
                    switch inp.MLI.method
                        case 'posneg'
                            predInterp = cell(nTs, nclass, 2);
                        case 'median'
                            predInterp = cell(nTs, nclass);
                    end
                    mapInterp = cell(nclass, nM);
                    mapInterp_ciu = cell(nclass, nM);
                    mapInterp_cil = cell(nclass, nM);
                    mapInterp_std = cell(nclass, nM);
                    for h=1:nclass
                        for nx = 1:nM
                            nY = size(inp.X(nx).Y,2);
                            mapInterp{h,nx} = zeros(nTs, nY);
                            mapInterp_ciu{h,nx} = zeros(nTs, nY);
                            mapInterp_cil{h,nx} = zeros(nTs, nY);
                            mapInterp_std{h,nx} = zeros(nTs, nY);
                        end
                    end

                    for q=1:numel(tInd) % Loop through CV2/OOCV cases 

                        fprintf('\n\n--- Working on case %g ---', tInd(q));
                        for nx = 1:numel(inp.X)
                            % Create artificial data and add it as "OOCV
                            % data" (inp.X(nx).Yocv) to the input structure. 
                            % if OOCV data has been selected for interpretation, then create
                            % inp.X(nx).Yocv2 and deal with this scenario in nk_ApplyTrainedPreproc 
                            % (see there)
                            Tr = inp.X(nx).Y(TrInd,:); covs = [];
                            if inp.oocvflag
                                if ~isempty(inp.covars_oocv), covs = inp.covars_oocv(tInd(q,:)); end
                                Ts = inp.X(nx).Yocv(tInd(q,:),:);
                            else
                                if ~isempty(inp.covars), covs = inp.covars(tInd(q,:)); end
                                Ts = inp.X(nx).Y(tInd(q,:),:);
                            end
                            inp = nk_CreateData4MLInterpreter( RandFeats(h, nx).I, Tr, Ts , covs, inp, nx );
                        end
                        
                        %% Step 4: generate predictions for artificial cases
                        for h=1:nclass  % Loop through binary comparisons

                            if nclass > 1, fprintf('\n*** %s #%g ***',algostr, h); end
    
                            [~, ~, nP, ~] = nk_GetModelParams2(analysis, multiflag, ll, h, inp.curlabel);
                           
                            switch inp.MLI.method
                                case 'posneg'
                                    inp.desc_oocv{1} = sprintf('%g%%-percentile modification', inp.MLI.upper_thresh);
                                    inp.desc_oocv{2} = sprintf('%g%%-percentile modification', inp.MLI.lower_thresh);
                                    [ inp, ~, ~, ~, ~, ~, ~, ~, mapYocv ] = nk_ApplyTrainedPreproc(analysis, inp, paramfl, Param);
                                case 'median'
                                    inp.desc_oocv = 'median modification';
                                    [ inp, ~, ~, ~, ~, ~, ~, ~, mapYocv] = nk_ApplyTrainedPreproc(analysis, inp, paramfl, Param);
                            end
                            for m = 1 : nP      % Loop through parameter combinations
                                for k=1:iy      % Loop through CV1 permutations
                                    for l=1:jy  % Loop through CV1 folds
                                    
                                        % Get feature feature subspace mask for current 
                                        % parameter grid position
                                        Fkl = GD.FEAT{Pspos(m)}{k,l,h};
    
                                        % Determine number of features in mask and
                                        % convert feature mask to logical index, if needed
                                        ul=size(Fkl,2); totLearn = totLearn + ul;
                                        if ~islogical(Fkl), F = Fkl ~= 0; else, F = Fkl; end

                                         % Set the pointer to the correct mapY shelf
                                        for n=1:numel(paramfl)
                                            pnt = 1;
                                            if ~BINMOD
                                                 if isfield(paramfl{n},'PREPROC') && ...
                                                   isfield(paramfl{n},'PXfull') && ...
                                                   ~isempty(paramfl{n}.P{1})
                                                    pnt = m;
                                                    break   
                                                end
                                            else
                                                if isfield(paramfl{n},'PREPROC') && ...
                                                   isfield(paramfl{n},'PXfull') && ...
                                                   ~isempty(paramfl{n}.P{h})
                                                    %Actual parameter node
                                                    pnt = m; 
                                                    break   
                                                end
                                            end
                                        end
  
                                        if BINMOD, hix = h; else, hix = 1; end
                                        [ TR , ~, ~, OCV ] = nk_ReturnAtOptPos(mapY.Tr{k,l}{hix},  mapY.CV{k,l}{hix}, mapY.Ts{k,l}{hix}, mapYocv.Ts{k,l}{hix}, Param{1}(k,l,hix), pnt);    

                                        switch inp.MLI.method

                                            case 'posneg'
                                                uD_pos = zeros(size(OCV{1},1),ul);
                                                uD_neg = zeros(size(OCV{1},1),ul);
                                                % Loop through feature subspaces
                                                for u=1:ul
            
                                                    % Extract features according to mask
                                                    TR_star   = nk_ExtractFeatures(TR, F, [], u);

                                                    % Apply trained model to
                                                    % artificial data and generate
                                                    % predictions for later
                                                    % evaluation
                                                    fprintf('\nCV2 [%g,%g], CV1 [%g,%g]: Model #%g => Predicting %g modified instances of case %g',  f, d, k, l, u, size(OCV{1},1), tInd(q));
                                                    [~, ~, uD_pos(:,u)] = nk_GetTestPerf(OCV{1}, ones(size(OCV{1},1),1), F(:,u), MD{h}{m}{k,l}{u}, TR_star, 1);
                                                    [~, ~, uD_neg(:,u)] = nk_GetTestPerf(OCV{2}, ones(size(OCV{2},1),1), F(:,u), MD{h}{m}{k,l}{u}, TR_star, 1);
        
                                                    % Detrend regressor predictions, if
                                                    % required by user input
                                                    if detrendfl 
                                                        beta = GD.Detrend{Pspos(m)}.beta;
                                                        p = GD.Detrend{Pspos(m)}.p;
                                                        uD_pos(:,u) = nk_DetrendPredictions2(beta, p, uD_pos(:,u));
                                                        uD_neg(:,u) = nk_DetrendPredictions2(beta, p, uD_neg(:,u));
                                                    end    
                                                end
                                                predInterp{q,h,1} = [predInterp{q,h,1} uD_pos];
                                                predInterp{q,h,2} = [predInterp{q,h,2} uD_neg];

                                            case 'median'

                                                uD = zeros(size(OCV,1),ul);
                                                % Loop through feature subspaces
                                                for u=1:ul
            
                                                    % Extract features according to mask
                                                    TR_star   = nk_ExtractFeatures(TR, F, [], u);

                                                    % Apply trained model to
                                                    % artificial data and generate
                                                    % predictions for later
                                                    % evaluation
                                                    fprintf('\nCV2 [%g,%g], CV1 [%g,%g]: Model #%g => Predicting %g modified instances of case %g', f, d, k, l, u, size(OCV,1), tInd(q));
                                                    [~, ~, uD(:,u)] = nk_GetTestPerf(OCV, ones(size(OCV,1),1), F(:,u), MD{h}{m}{k,l}{u}, TR_star, 1);
        
                                                    % Detrend regressor predictions, if
                                                    % required by user input
                                                    if detrendfl 
                                                        beta = GD.Detrend{Pspos(m)}.beta;
                                                        p = GD.Detrend{Pspos(m)}.p;
                                                        uD(:,u) = nk_DetrendPredictions2(beta, p, uD(:,u)); 
                                                    end    
                                                end
                                                predInterp{q,h} = [predInterp{q,h} uD];
                                        end  
                                    end
                                end
                            end
                        end
                       
                        % Remove artificial data from inp structure
                        % otherwise nk_ApplyTrainedPreproc may not work
                        % properly.
                        if inp.oocvflag
                            inp.X = rmfield(inp.X,'Yocv2');
                        else
                            inp.X = rmfield(inp.X,'Yocv');
                        end
                        
                        %% Step 5: Evaluate impact of input data modifications using obtained predictions
                        for h=1:nclass
                            Oh = nm_nanmedian(predOrig{h}(q,:),2);
                            for nx = 1:nM
                                switch inp.MLI.method
                                    case 'posneg'
                                        Rh = [ nm_nanmedian(predInterp{q,h,1},2) nm_nanmedian(predInterp{q,h,2},2)]; 
                                    case 'median'
                                        Rh = nm_nanmedian(predInterp{q,h},2);
                                end
                            end
                            [mapInterp{h, nx}(q,:), ...
                                mapInterp_ciu{h, nx}(q,:), ...
                                mapInterp_cil{h, nx}(q,:), ...
                                mapInterp_std{h, nx}(q,:)] = nk_MapModelPredictions(Oh, Rh, inp.X(nx).I, inp.MLI.method, inp.MLI.RangePred(h));
                        end
                    end
                    fprintf('\nSaving %s', oMLIpath); 
                    save(oMLIpath,'predOrig', 'predInterp', 'mapInterp', 'mapInterp_ciu', 'mapInterp_cil', 'mapInterp_std', 'operm','ofold');
                    if saveparam 
                        fprintf('\nSaving %s', OptModelPath); save(OptModelPath, 'MD', 'ofold','operm'); 
                    end
                end
                
            case 1
                
                vpth = deblank(matfiles{f,d});
                if isempty(vpth) || ~exist(vpth,'file') && GridAct(f,d)
                    error(['No valid MLIdatamat detected for CV2 partition ' '[' num2str(f) ', ' num2str(d) ']!']);
                else
                    [~,vnam] = fileparts(vpth);
                    fprintf('\n\nLoading MLI results for CV2 partition [ %g, %g ]:', f, d);
                    fprintf('\n%s',vnam);
                    load(vpth)
                end 
        end
        
        %% Step 6: Concatenate results and assign them to results container
        [RootPath, FileNames{f,d}] = fileparts(oMLIpath);  
        ll=ll+1;
        for h = 1:nclass
            for nx = 1 : nM
                switch MODEFL
                    case 'classification'
                        Results.BinResults(h).Modality(nx).Y_mapped(tInd,:,f) = mapInterp{h,nx};
                        Results.BinResults(h).Modality(nx).Y_mapped_ciu(tInd,:,f) = mapInterp_ciu{h,nx};
                        Results.BinResults(h).Modality(nx).Y_mapped_cil(tInd,:,f) = mapInterp_cil{h,nx};
                        Results.BinResults(h).Modality(nx).Y_mapped_std(tInd,:,f) = mapInterp_std{h,nx};
                    case 'regression'
                        Results.RegrResults.Modality(nx).Y_mapped(tInd,:,f) = mapInterp{h,nx}/ix;
                        Results.RegrResults.Modality(nx).Y_mapped_ciu(tInd,:,f) = mapInterp_ciu{h,nx};
                        Results.RegrResults.Modality(nx).Y_mapped_cil(tInd,:,f) = mapInterp_cil{h,nx};
                        Results.RegrResults.Modality(nx).Y_mapped_std(tInd,:,f) = mapInterp_std{h,nx};
                end
            end
        end
    end
    ol=ol+1;
end
ol=ol-1;
for h = 1:nclass
    for nx = 1 : nM
        switch MODEFL
            case 'classification'
                Results.BinResults(h).Modality(nx).Y_mapped     = nm_nanmean(Results.BinResults(h).Modality(nx).Y_mapped(:,:,1:ol),3);
                Results.BinResults(h).Modality(nx).Y_mapped_ciu = nm_nanmean(Results.BinResults(h).Modality(nx).Y_mapped_ciu(:,:,1:ol),3);
                Results.BinResults(h).Modality(nx).Y_mapped_cil = nm_nanmean(Results.BinResults(h).Modality(nx).Y_mapped_cil(:,:,1:ol),3);
                Results.BinResults(h).Modality(nx).Y_mapped_std = nm_nanmean(Results.BinResults(h).Modality(nx).Y_mapped_std(:,:,1:ol),3);
            case 'regression'
                Results.RegrResults.Modality(nx).Y_mapped       = nm_nanmean(Results.RegrResults.Modality(nx).Y_mapped(:,:,1:ol),3);
                Results.RegrResults.Modality(nx).Y_mapped_ciu   = nm_nanmean(Results.RegrResults.Modality(nx).Y_mapped_ciu(:,:,1:ol),3);
                Results.RegrResults.Modality(nx).Y_mapped_cil   = nm_nanmean(Results.RegrResults.Modality(nx).Y_mapped_cil(:,:,1:ol),3);
                Results.RegrResults.Modality(nx).Y_mapped_std   = nm_nanmean(Results.RegrResults.Modality(nx).Y_mapped_std(:,:,1:ol),3);
        end
    end
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

