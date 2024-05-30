function [Results, FileNames, RootPath] = nk_OOCV(inp)
% =========================================================================
% FORMAT Results = nk_OOCV(inp)
% =========================================================================
% Independent test data prediction module
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, last modified 11/2023

global SVM RFE MULTI MODEFL CV EVALFUNC OOCV SCALE SAV CVPOS RAND
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
FullPartFlag    = RFE.ClassRetrain;
switch inp.analmode
    case 0
        ovrwrt  = inp.ovrwrt;                       % overwrite existing data files
    case 1
        oocvmat  = inp.oocvmat;                     % OOCVdatamat
end
multiflag       = false; if strcmp(MODEFL,'classification'), multiflag = inp.multiflag; end
TrainWithCV2Ts  = OOCV.trainwithCV2Ts;
saveparam       = inp.saveparam;
loadparam       = inp.loadparam;
nclass          = inp.nclass;
ngroups         = inp.ngroups;
nsubgroups      = inp.nsubgroups;
analysis        = inp.analysis;
GridAct         = inp.GridAct;
batchflag       = inp.batchflag;
algostr         = GetMLType(SVM);
[ylm, ylb]      = nk_GetScaleYAxisLabel(SVM);

% Setup CV2 and CV1 counters and data containers:
[ix, jx]        = size(CV.TrainInd);
[iy, jy]        = size(CV.cvin{1,1}.TrainInd);
binOOCVD        = cell(nclass,1);
ll              = 1; 
totLearn        = 0;
ha              = [];
if ~exist('GridAct','var') || isempty(GridAct), GridAct = nk_CVGridSelector(ix,jx); end
if ~exist('batchflag','var') || isempty(batchflag), batchflag = false; end

CLtarg_m = '*-';
CLdec_m = 'o-';
CL = getNMcolors;

% Train models with CV2 test data
switch TrainWithCV2Ts
    case 1
        if FullPartFlag
            % Use entire inner-cycle data to retrain models
            TrainWithCV2TsStr = 'CV1-Tr + CV1-Ts';
        else
            % Use only inner-cycle training to retrain models
            TrainWithCV2TsStr = 'CV1-Tr';
        end
    case 2
        % Use all available data t
        TrainWithCV2TsStr = 'CV1-Tr + CV1-Ts + CV2-Ts';
end

% Check and transform labels if needed
% Also if we are in multilabel mode this function return the current label
% in inp.label
inp = nk_ApplyLabelTransform( SCALE, MODEFL, inp );

% Check whether oocv labels are available and process them accordingly.
% Check whether given analysis used alternative label and whether it's
% available
if ~isfield(inp,'labelOOCV')
    labelOOCV = zeros(inp.nOOCVsubj,1);
    LabelMode = false;
else
    labelOOCV   = inp.curlabelOOCV;
    LabelMode   = true;
    if ~isempty(inp.PolyFact), labelOOCV = labelOOCV .^ (1/inp.PolyFact); end
    indDicho    = cell(nclass,1);
    labelDicho  = cell(nclass,1);
    nCV2 = sum(GridAct(:));
    if strcmp(MODEFL,'classification')
        Results.BinCV2Performance_DecisionValues_History    = nan(nclass,nCV2 );
        Results.BinCV2Performance_Targets_History           = nan(nclass,nCV2 );
        Results.BinCV2Performance_DecisionValues            = nan(nclass,nCV2 );
        Results.BinCV2Performance_Targets                   = nan(nclass,nCV2 );
        % Create binary indices
        for curclass=1:nclass
            if numel(CV.class{1,1}{curclass}.groups)>1
                indClass1 = (labelOOCV == CV.class{1,1}{curclass}.groups(1));
                indClass2 = (labelOOCV == CV.class{1,1}{curclass}.groups(2));
            else
                indClass1 = labelOOCV == curclass; indClass2 = ~indClass1;
            end
            indDicho{curclass}              = indClass1 | indClass2;
            labelDicho{curclass}            = zeros(numel(indDicho{curclass}),1);
            labelDicho{curclass}(indClass1) = 1; labelDicho{curclass}(indClass2) = -1;
        end
        if inp.PERM.flag==1
            Results.PermAnal.ModelPermPerf = cell(nclass,1);
            Results.PermAnal.ModelObsPerf(curclass) = zeros(nclass,1);
            Results.PermAnal.ModelPermPerfCrit{curclass} = cell(nclass,1);
            Results.PermAnal.ModelPermSignificance(curclass) = zeros(nclass,1); 
            compfun = str2func(nk_ReturnEvalOperator(SVM.GridParam));
        end
    else
        Results.CV2Performance_PredictedValues_History = nan(1, nCV2 );
        Results.CV2Performance_PredictedValues = nan(1, nCV2 );
        labelDicho{1} = ones(size(inp.labelOOCV));
        if inp.PERM.flag==1
             compfun = str2func(nk_ReturnEvalOperator(SVM.GridParam));
        end
    end
end

if MULTI.flag && multiflag
    Results.MultiCV2Performance_History = [];
    Results.MultiCV2PerformanceLL = [];
    Results.MultiCV2PredictionsLL = [];
end

% Check whether you have to perform permutation analysis on OOCV data and
% create permutation index file
if inp.PERM.flag==1
    if isfield(inp.PERM,'permdeffile') && ~isempty(inp.PERM.permdeffile) && exist(inp.PERM.permdeffile,'file')
        fprintf('\nLoading permutation definition: %s', inp.PERM.permdeffile);
        load(inp.PERM.permdeffile)
        if ~exist("indpermA","var")
            error('\n%s does not contain permutation definitions!', inp.PERM.permdeffile)
        end
        inp.PERM.nperms = size(indpermA,2);
        fprintf('\nn=%g permutations found in permutation matrix.',inp.PERM.nperms);
    else
        inp.PERM.suffix = '_OOCVpermmat_ID';
        indpermA = nk_GenPermMatrix(CV, SAV, inp); 
    end
    binOOCVD_perm = cell(nclass,1);
    for h=1:nclass
        if ~RFE.CV2Class.EnsembleStrategy.AggregationLevel
            binOOCVD_perm{h} = nan(inp.nOOCVsubj, sum(inp.analysis.NumModels(h,:, inp.curlabel)>0), inp.PERM.nperms);
        else
            binOOCVD_perm{h} = nan(inp.nOOCVsubj, sum(inp.analysis.NumModels(h, :, inp.curlabel)), inp.PERM.nperms); % this can cause oom issues if many observations, models and permutaitons
        end
    end
end

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
if isfield(RAND,'Decompose') && RAND.Decompose == 2
    BINMOD = 0;
end

CVPOS.fFull = FullPartFlag;
FileNames = cell(ix,jx); fnd=false;

% Parameter flag structure for preprocessing
paramfl = struct('use_exist',inp.loadparam, ...
                 'found', false, ...
                 'write', inp.saveparam, ...
                 'writeCV1', inp.saveCV1, ...
                 'multiflag', multiflag);

%Pre-smooth data, if needed, to save computational time
inp.ll=inp.GridAct';inp.ll=find(inp.ll(:));
if ~inp.analmode
    inp = nk_PerfInitSpatial(analysis, inp, paramfl);
end
cntOOCVDh = zeros( ix*jx, nclass, 2);
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
        CVPOS.CV2p = f;
        CVPOS.CV2f = d;
        
        binOOCVDh = cell(nclass,1); inp.f = f; inp.d = d; inp.ll = ll;
        if inp.PERM.flag == 1
            binOOCVDh_perm = cell(nclass,1); 
            for h=1:nclass
                binOOCVDh_perm{h} = zeros(inp.nOOCVsubj, inp.analysis.NumModels(h, ll, inp.curlabel), inp.PERM.nperms);
            end
        end
        operm = f; ofold = d;
        % Create OOCV partition file path
        oOOCVpath = nk_GenerateNMFilePath(inp.rootdir, SAV.matname, inp.datatype, inp.multlabelstr, inp.varstr, inp.id, operm, ofold);
        OptModelPath = nk_GenerateNMFilePath( inp.saveoptdir, SAV.matname, 'OptModel', [], inp.varstr, inp.id, operm, ofold);
        switch inp.analmode
            case 0
                 loadfl = false;
                 if exist(oOOCVpath,'file') && ~ovrwrt && ~batchflag
                    
                    [~, onam] = fileparts(oOOCVpath);
                    fprintf('\nOOCVdatamat found for CV2 [%g,%g]:',f,d)
                    fprintf('\nLoading: %s',onam)
                    try
                        load(oOOCVpath)
                        loadfl = true;
                    catch
                        fprintf('\nCould not open file. May be corrupt. Recompute CV2 partition [%g,%g].',f,d);
                        loadfl = false;
                    end
                    
                elseif exist(oOOCVpath,'file') && batchflag
                     % in batch mode we do not compute statistics across the
                    % CV2 partitions
                    [~, onam] = fileparts(oOOCVpath);
                    fprintf('\nOOCVdatamat found for CV2 [%g,%g]:\n%s',f,d,onam)
                    fprintf('\nBatch mode detected. Continue.')
                    [RootPath, FileNames{f,d}] = fileparts(oOOCVpath); 
                    load(oOOCVpath)
                    loadfl=true;
                end
                   
                if ~loadfl 

                    [ inp, contfl, analysis, mapY, GD, MD, Param, paramfl, mapYocv ] = nk_ApplyTrainedPreproc(analysis, inp, paramfl);

                    if contfl, continue; end

                    fndMD = false; 
                    if loadparam && isfield(inp,'optmodelmat') && exist(inp.optmodelmat{operm,ofold},'file')
                        fprintf('\nLoading OptModel: %s', inp.optmodelmat{operm,ofold});
                        load(inp.optmodelmat{operm,ofold},'MD'); fndMD = true; 
                    end
                    if ~fndMD, MD = cell(nclass,1); end

                    % -----------------------------------------------------------------
                    for h=1:nclass % Loop through binary comparisons

                        if nclass > 1, fprintf('\n\n*** %s #%g ***',algostr, h); end

                        %% Step 1: Get optimal model parameters
                        % Retrieve optimal parameters from precomputed analysis structure
                        % Differentiate according to binary or multi-group mode
                        [Ps, Pspos, nP, Pdesc] = nk_GetModelParams2(analysis, multiflag, ll, h, inp.curlabel);

                        if ~fndMD , MD{h} = cell(nP,1); end
                        cntModels = 1;

                        % Loop through parameter combinations
                        for m = 1 : nP

                            cPs = Ps(m,:); sPs = nk_PrepMLParams(Ps, Pdesc, m);
                            P_str = nk_DefineMLParamStr(cPs, analysis.Model.ParamDesc, h);

                            %% Step 2: Apply trained model to OOCV data 
                            % Optionally, retrain every base learner in current CV1
                            % [k,l] partition (across different grid positions, if available)
                            if ~fndMD,MD{h}{m} = cell(iy,jy); end

                            for k=1:iy % CV1 permutations

                                for l=1:jy % CV1 folds
                                    
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
                                    CVInd = mapY.CVInd{k,l}{h};
                                    TrInd = mapY.TrInd{k,l}{h};
                                    TsInd = mapY.TsInd{h};

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

                                    % get training data using pointers
                                    % Either (a) only CV1 training data, (b) CV1
                                    % training and test data, (c) CV1 training & test
                                    % data as well as CV2 test data              
                                    if BINMOD, hix = h; else, hix = 1; end
                                    [ TR , CV1, ~, OCV ] = nk_ReturnAtOptPos(mapY.Tr{k,l}{hix},  mapY.CV{k,l}{hix}, mapY.Ts{k,l}{hix}, mapYocv.Ts{k,l}{hix}, Param{1}(k,l,hix), pnt);                                        
                                    modelTrL = mapY.TrL{k,l}{h}; 

                                    if FullPartFlag
                                        % Merge CV1 training and CV1 test data if recomputing 
                                        % of model using entire CV1 partition data had been previously used
                                        TR = [ TR; CV1 ]; 
                                        % Merge CV1 training and CV1 test data
                                        modelTrL = [modelTrL; mapY.CVL{k,l}{h}]; 
                                        % Merge CV1 training and test indices
                                        TrInd = [TrInd; CVInd]; 
                                    end
                                  
                                    % Prepare decision / probability value container
                                    uD = zeros(inp.nOOCVsubj,ul);

                                    % Loop through feature subspaces
                                    if ~fndMD 
                                        MD{h}{m}{k,l} = cell(ul,1); 
                                        fprintf(['\nRetrain models in CV2 [%g,%g], ' ...
                                            'CV1 [%g,%g], %g %s, (total # learners: %g) => Data: %s, ML params [%s] ==> '], ...
                                            f, d, k, l, ul, algostr, totLearn, TrainWithCV2TsStr, P_str)

                                        % Impute labels if needed
                                        [modelTrL, TR, TrInd] = nk_LabelImputer( modelTrL, TR, TrInd, sPs, IMPUTE);
                                        TR = TR(TrInd,:);
                                    else
                                        fprintf(['\nUse precomputed models in CV2 [%g,%g], ' ...
                                            'CV1 [%g,%g], %g %s, (total # learners: %g) => Data: %s, ML params [%s] ==> '], ...
                                            f, d, k, l, ul, algostr, totLearn, TrainWithCV2TsStr, P_str)
                                        fprintf('Apply to OOCV data. ');
                                    end

                                    if inp.PERM.flag==1
                                        if iscell(indpermA)
                                            indperm = indpermA{h}(TrInd,:); 
                                        else
                                            indperm = indpermA(TrInd,:); 
                                        end
                                        uD_perm = zeros(inp.nOOCVsubj,ul,inp.PERM.nperms);
                                    end

                                    % Loop through feature subspaces
                                    for u=1:ul

                                        % Extract features according to mask
                                        Ymodel = nk_ExtractFeatures(TR, F, [], u);
                                        
                                        if ~fndMD
                                            fprintf('Computing OptModel');
                                            [~, MD{h}{m}{k,l}{u}] = nk_GetParam2(Ymodel, modelTrL, sPs, 1);
                                        end

                                        % Apply model to independent test data 
                                        [~, ~, uD(:,u)] = nk_GetTestPerf(OCV, labelOOCV, F(:,u), MD{h}{m}{k,l}{u}, Ymodel, 1);
                                       
                                        % Detrend regressor predictions, if
                                        % required by user input
                                        if detrendfl
                                            beta = GD.Detrend{Pspos(m)}.beta;
                                            p = GD.Detrend{Pspos(m)}.p;
                                            uD(:,u) = nk_DetrendPredictions2(beta, p, uD(:,u)); 
                                        end

                                        if inp.PERM.flag==1
                                            fprintf('\nRunning permutation analysis\t');
                                            % retrieve permuted labels
                                            % retrain model using permuted
                                            % labels
                                            for curperm = 1:inp.PERM.nperms
                                                fprintf('.')
                                                modelTrL_perm = nk_VisXPermY(Ymodel, inp.label, TrInd, 1, indperm, [], curperm);
                                                [ ~, MD_perm ] = nk_GetParam2(Ymodel, modelTrL_perm , sPs, 1);
                                                [ ~, ~, uD_perm( :, u, curperm) ] = nk_GetTestPerf(OCV, labelOOCV, F(:,u), MD_perm, Ymodel, 1);
                                                if detrendfl
                                                    uD_perm(:,u, curperm) = nk_DetrendPredictions2(beta, p, uD_perm(:,u, curperm)); 
                                                end
                                            end
                                        end
                                    end

                                    fprintf('. Done');
                                    
                                    %% Step 3: Concatenate binary classifiers / predictors into [k,l,h] array
                                    binOOCVDh{h} = [binOOCVDh{h} uD];
                                    if inp.PERM.flag == 1
                                       binOOCVDh_perm{h}(:,cntModels:cntModels+size(uD,2)-1,:) = uD_perm;
                                    end
                                    cntModels = cntModels + size(uD,2);
                                end
                            end
                        end
                        % Bug fix for wrong indexing causing nk_OOCV to
                        % crash in permutation mode when concatenating
                        % all predictions instead of concatenating the
                        % CV1 mean predictions.
                        if ll==1
                            cntOOCVDh(ll, h, 1) = 1;
                            cntOOCVDh(ll, h, 2) = cntModels-1;
                        else
                            cntOOCVDh(ll, h, 1) = cntOOCVDh(ll-1, h, 2) + 1;
                            cntOOCVDh(ll, h, 2) = cntOOCVDh(ll, h, 1) + cntModels-2;
                        end
                    end
                    fprintf('\nSaving %s', oOOCVpath); 
                    if inp.PERM.flag==1
                       save(oOOCVpath,'binOOCVDh','binOOCVDh_perm','cntOOCVDh','operm','ofold');
                    else
                       save(oOOCVpath,'binOOCVDh','cntOOCVDh','operm','ofold');
                    end
                    if saveparam, fprintf('\nSaving %s', OptModelPath); save(OptModelPath, 'MD', 'ofold','operm'); end
                end
                
            case 1
                
                vpth = deblank(oocvmat{f,d});
                if isempty(vpth) || ~exist(vpth,'file') && GridAct(f,d)
                    error(['No valid OOCVdatamat detected for CV2 partition ' '[' num2str(f) ', ' num2str(d) ']!']);
                else
                    [~,vnam] = fileparts(vpth);
                    fprintf('\n\nLoading independent test data results for CV2 partition [ %g, %g ]:', f, d);
                    fprintf('\n%s',vnam);
                    load(vpth)
                end 

        end
        
        for curclass=1:nclass
            [~, ~, nP, ~] = nk_GetModelParams2(analysis, multiflag, ll, curclass, inp.curlabel);
            if ~RFE.CV2Class.EnsembleStrategy.AggregationLevel
                fprintf('\nCompute mean of base learners'' outputs of current CV2 partition and add them to the ensemble matrix.')
                EnsDat = nm_nanmedian(binOOCVDh{curclass},2);
                if inp.PERM.flag==1
                    for curperm=1:inp.PERM.nperms
                        binOOCVD_perm{curclass}(:,ll,curperm) = nm_nanmedian(binOOCVDh_perm{curclass}(:,:,curperm),2);
                    end
                end
            else
                EnsDat = binOOCVDh{curclass};
                fprintf('\nAdd all base learners'' outputs to the ensemble matrix without averaging.')
                if inp.PERM.flag == 1 
                    for m=1:nP
                        for curperm=1:inp.PERM.nperms
                            binOOCVD_perm{curclass}(:,cntOOCVDh(ll, curclass, 1):cntOOCVDh(ll, h, 2),curperm) = ...
                                binOOCVDh_perm{curclass}(:,:,curperm);
                        end
                    end
                end
            end
            binOOCVD{curclass} = [binOOCVD{curclass} EnsDat];
        end
        
        %% Step 4: Compute OOCV multi-group prediction from current binary classifier arrays
        if MULTI.flag && multiflag

            %% Step4a: Get multi-group prediction for current CV2 partition
            multiOOCVll = []; multiClassll = [];
            for curclass = 1 : nclass
                multiOOCVll = [multiOOCVll binOOCVDh{curclass}];
                multiClassll = [multiClassll ones(1,size(binOOCVDh{curclass},2))*curclass];
            end

            % Compute multi-group labels (& performance, if labels are
            % available)
            [MultiCV2PerformanceLL, MultiCV2PredictionsLL] = ...
                nk_MultiEnsPerf(multiOOCVll, sign(multiOOCVll), labelOOCV, multiClassll, ngroups);

            if LabelMode
                Results.MultiCV2PerformanceLL = [ Results.MultiCV2PerformanceLL MultiCV2PerformanceLL ]; 
            end
            Results.MultiCV2PredictionsLL = [ Results.MultiCV2PredictionsLL MultiCV2PredictionsLL ];
            
            %% Step 4b: Compute multi-group prediction based on ensemble generated from ...
            % all base learners (from the start to the current position of
            % the CV2 loop)
            nDicho = 0; nDichoH = zeros(nclass+1,1); nDichoH(1) = 1;
            for curclass = 1: nclass
                nDicho = nDicho + size(binOOCVD{curclass},2);
                nDichoH(curclass+1) = 1 + nDicho;
            end
            multiOOCV  = zeros(inp.nOOCVsubj, nDicho);
            multiClass = zeros(1, nDicho);
            for curclass=1:nclass
                multiOOCV( : , nDichoH(curclass) : nDichoH(curclass+1) - 1) = binOOCVD{curclass};
                multiClass(1, nDichoH(curclass) : nDichoH(curclass+1) - 1) = ones(1,size(binOOCVD{curclass},2)) * curclass;
            end

            % Compute multi-group labels (& performance, if labels are
            % available)
            [Results.MultiCV2Performance, Results.MultiCV2Predictions] = ...
                nk_MultiEnsPerf(multiOOCV, sign(multiOOCV), labelOOCV, multiClass, ngroups);
            if LabelMode
                Results.MultiCV2Performance_History = ...
                    [ Results.MultiCV2Performance_History Results.MultiCV2Performance];
            end

        end

        indnan = isnan(labelOOCV) | sum(isnan(binOOCVD{1}),2)==size(binOOCVD{1},2);

        if LabelMode
        
            %% Step 5: Assess binary classifier performance, if OOCV Label has been specified
            if strcmp(MODEFL,'classification')
                for curclass=1:nclass
                    
                    % These are the algorithm outputs at the current CV2
                    % partition which are stored in binOOCVDh
                    binOOCVDhx_ll = binOOCVDh{curclass}(indDicho{curclass},:);
                    
                    % ... and binOOCVD contains the outputs produced by the
                    % ensemble across the CV2 partitions processed so far
                    binOOCVDhx = binOOCVD{curclass}(indDicho{curclass},:);

                    if size(binOOCVDhx,2)>1
                        hrx = sign(nm_nansum(sign(binOOCVDhx),2)); 
                        hdx = nm_nanmedian(binOOCVDhx,2); hdx(indnan) = nan;
                    else
                        hrx = sign(binOOCVDhx);
                        hdx = binOOCVDhx;
                    end
                    if size(binOOCVDhx_ll,2)>1
                        hrx_ll = sign(nm_nansum(sign(binOOCVDhx_ll),2)); 
                        hdx_ll = nm_nanmedian(binOOCVDhx_ll,2); 
                    else
                        hrx_ll = sign(binOOCVDhx_ll);
                        hdx_ll = binOOCVDhx_ll;
                    end
                    hdx(indnan) = nan;
                    hdx_ll(indnan) = nan;
                    if ll==1 
                        Results.BinLabels{curclass} = labelDicho{curclass}; 
                        if inp.multiflag, Results.Labels = labelOOCV; end
                    end
                    Results.BinCV2Performance_Targets(curclass, ll) = EVALFUNC(labelDicho{curclass}(indDicho{curclass}) , hrx_ll);
                    Results.BinCV2Performance_DecisionValues(curclass,ll) = EVALFUNC(labelDicho{curclass}(indDicho{curclass}), hdx_ll);
                    Results.BinCV2Sensitivity_DecisionValues(curclass,ll) = SENSITIVITY(labelDicho{curclass}(indDicho{curclass}), hdx_ll);
                    Results.BinCV2Specificity_DecisionValues(curclass,ll) = SPECIFICITY(labelDicho{curclass}(indDicho{curclass}), hdx_ll);
                    Results.contingency{curclass} = ALLPARAM(labelDicho{curclass}(indDicho{curclass}), hdx);
                    Results.BinCV2Performance_Targets_History(curclass,ll) = EVALFUNC(labelDicho{curclass}(indDicho{curclass}) , hrx);
                    Results.BinCV2Performance_DecisionValues_History(curclass,ll) = EVALFUNC(labelDicho{curclass}(indDicho{curclass}) , hdx);

                    % Subgroup information available?
                    % If so, process in one of two modes:
                    % mode 1: no reference group provided, performance
                    % metrics are computed for each subgroup
                    % mode 2: reference group provided, performance metrics
                    % are computed for each group against reference group
                    % (e.g. HCs subjects)
                    if isfield(inp,'groupind')
                        for g = 1:nsubgroups
                            indr = []; 
                            if iscell(inp.groupvec)
                                indg = find(strcmp(inp.groupind,inp.groupvec{g}));
                            else
                                indg = find(inp.groupind == inp.groupvec(g));
                            end
                            refstr = sprintf('%s', inp.groupnames{g});
                            if isfield(inp,'refgroup')
                                if iscell(inp.groupvec)
                                    indr = find(strcmp(inp.groupind,inp.refgroup));
                                else
                                    indr = find(inp.groupind == inp.refgroup);
                                end
                                refstr = sprintf('%s vs %s', inp.groupnames{g}, inp.groupnames{inp.refgroup});
                            end
                            idx = [indg; indr];
                            Results.Group{g}.GroupName = refstr; 
                            Results.Group{g}.ObservedValues{curclass} = labelDicho{curclass}(idx);
                            Results.Group{g}.BinCV2Performance_Targets(curclass, ll) = EVALFUNC(Results.Group{g}.ObservedValues{curclass}(indDicho{curclass}(idx)) , hrx_ll(idx));
                            Results.Group{g}.BinCV2Performance_DecisionValues(curclass, ll) = EVALFUNC(Results.Group{g}.ObservedValues{curclass}(indDicho{curclass}(idx)), hdx_ll(idx));
                            Results.Group{g}.BinCV2Performance_Targets_History(curclass,ll) = EVALFUNC(Results.Group{g}.ObservedValues{curclass}(indDicho{curclass}(idx)) , hrx(idx));
                            Results.Group{g}.BinCV2Performance_DecisionValues_History(curclass,ll) = EVALFUNC(Results.Group{g}.ObservedValues{curclass}(indDicho{curclass}(idx)) , hdx(idx));
                        end
                    end
                end
            else % regression case
                if size(binOOCVDh{1},2)>1
                    hdx_ll = nm_nanmedian(binOOCVDh{1},2); 
                else
                    hdx_ll = binOOCVDh{1}; 
                end
                if size(binOOCVD{1},2)>1
                    hdx = nm_nanmedian(binOOCVD{1},2); 
                else
                    hdx = binOOCVD{1}; 
                end
                hdx_ll(indnan) = nan;
                hdx(indnan) = nan;
                if ll==1 || batchflag, Results.RegrLabels = labelOOCV; Results.RegrLabels(indnan) = nan; end
                Results.CV2Performance_PredictedValues(ll) = EVALFUNC(Results.RegrLabels, hdx_ll);
                Results.CV2Performance_PredictedValues_History(ll) = EVALFUNC(Results.RegrLabels, hdx);
                if isfield(inp,'groupind')
                     for g = 1:nsubgroups
                        if iscell(inp.groupvec)
                            indg = find(strcmp(inp.groupind,inp.groupvec{g}));
                        else
                            indg = find(inp.groupind == inp.groupvec(g));
                        end
                        idx = indg;
                        Results.Group{g}.GroupName = inp.groupnames{g};
                        Results.Group{g}.CV2Performance_PredictedValues(ll) = EVALFUNC(Results.RegrLabels(idx), hdx_ll(idx));
                        Results.Group{g}.CV2Performance_PredictedValues_History(ll) = EVALFUNC(Results.RegrLabels(idx), hdx(idx));
                     end
                end
            end
            % Disply progress information in figure
            if ~batchflag
                if ~exist('hu','var') || isempty(hu) 
                    hu = findobj('Tag','OOCV');
                    if isempty(hu)
                        sz = get(0,'ScreenSize');
                        win_wdth = sz(3)/2; win_hght = sz(4)/1.25; win_x = sz(3)/2 - win_wdth/2; win_y = sz(4)/2 - win_hght/2;
                        hu = figure('Tag','OOCV', ...
                            'NumberTitle','off', ...
                             'MenuBar','none', ...
                             'Color', [0.9 0.9 0.9], ...
                             'Position', [win_x win_y win_wdth win_hght]); 
                    end
                end
                hu.Name = sprintf('NM Application Viewer: %s => %s [OOCV #%g]',inp.analysis_id, inp.desc_oocv, inp.oocvind);
                set(0,'CurrentFigure',hu); if ll==1, clf ; end
                if ll==1 || isempty(ha); ha = tight_subplot(2,nsubgroups,0.03,[.1 .05],[.1 .05]); end
                lg = cell(nsubgroups,nclass); hc = zeros(nsubgroups, nclass); lh = lg; hd = hc;
                for curgroup=1:nsubgroups
                    hold(ha(curgroup),'on');hold (ha(curgroup+nsubgroups),'on');
                    if ll==1 && nsubgroups>1
                        title(ha(curgroup), Results.Group{curgroup}.GroupName, 'FontWeight', 'bold', 'Interpreter', 'none');
                    end
                    for curclass=1:nclass
                        switch MODEFL
                            case 'classification'
                                if nsubgroups>1
                                    yvalTH = Results.Group{curgroup}.BinCV2Performance_Targets_History(curclass,:);
                                    yvalDH = Results.Group{curgroup}.BinCV2Performance_DecisionValues_History(curclass,:);
                                    yvalTC = Results.Group{curgroup}.BinCV2Performance_Targets(curclass,:);
                                    yvalDC = Results.Group{curgroup}.BinCV2Performance_DecisionValues(curclass,:);
                                    lgH = cellstr([sprintf("Cl#%g: target perf [overall]",curclass); ... 
                                                        sprintf("Cl#%g: score perf [overall]",curclass)]); 
                                    lgC = cellstr([sprintf("Cl#%g: target perf [current]",curclass); ...
                                                        sprintf("Cl#%g: score perf [current]",curclass)]); 
                                else
                                    yvalTH = Results.BinCV2Performance_Targets_History(curclass,:);
                                    yvalDH = Results.BinCV2Performance_DecisionValues_History(curclass,:);
                                    yvalTC = Results.BinCV2Performance_Targets(curclass,:);
                                    yvalDC = Results.BinCV2Performance_DecisionValues(curclass,:);
                                    lgH = cellstr([sprintf("Cl#%g: target perf [overall]",curclass); ... 
                                                        sprintf("Cl#%g: score perf [overall]",curclass)]); 
                                    lgC = cellstr([sprintf("Cl#%g: target perf [current]",curclass); ...
                                                        sprintf("Cl#%g: score perf [current]",curclass)]); 
                                end
                                hc(curgroup,curclass) = plot(ha(curgroup),yvalTH,CLtarg_m,'Color', CL(curclass,:));
                                hc(curgroup,curclass) = plot(ha(curgroup),yvalDH, CLdec_m,'Color', CL(curclass,:));
                                hd(curgroup,curclass) = plot(ha(curgroup+nsubgroups),yvalTC,CLtarg_m,'Color', CL(curclass,:));
                                hd(curgroup,curclass) = plot(ha(curgroup+nsubgroups),yvalDC, CLdec_m,'Color', CL(curclass,:));
                            
                                lg{curgroup} = [lg{curgroup}; lgH];
                                lh{curgroup} = [lh{curgroup}; lgC];
                                
                            case 'regression'
                                if nsubgroups>1
                                    hd(curgroup,1) = plot(ha(curgroup),Results.Group{curgroup}.CV2Performance_PredictedValues_History, CLtarg_m, 'Color', CL(curclass,:));
                                    lg{curgroup,1} = sprintf('Predictor performance of current ensemble [overall]'); 
                                    hd(curgroup,1) = plot(ha(curgroup+nsubgroups),Results.Group{curgroup}.CV2Performance_PredictedValues, CLtarg_m, 'Color', CL(curclass,:));
                                    lh{curgroup,1} = sprintf('Predictor performance of current ensemble [current]'); 
                                else
                                    hd(curgroup,1) = plot(ha(curgroup),Results.CV2Performance_PredictedValues_History, CLtarg_m, 'Color', CL(curclass,:));
                                    lg{curgroup,1} = sprintf('Predictor performance of current ensemble [overall]'); 
                                    hd(curgroup,1) = plot(ha(curgroup+nsubgroups),Results.CV2Performance_PredictedValues, CLtarg_m, 'Color', CL(curclass,:));
                                    lh{curgroup,1} = sprintf('Predictor performance of current ensemble [current]'); 
                                end
                        end
                    end
                    if MULTI.flag && multiflag,plot(Results.MultiCV2Performance_History,'k+'); end
                    ylim(ha(curgroup),ylm); ha(curgroup).YTickLabelMode='auto';
                    ylim(ha(curgroup+nsubgroups),ylm); ha(curgroup+nsubgroups).YTickLabelMode='auto';
                    if ll_start+nCV2-1 > ll_start
                        xlim(ha(curgroup),[ll_start ll_start+nCV2-1]); ha(curgroup).XTickLabelMode='auto';
                        xlim(ha(curgroup+nsubgroups),[ll_start ll_start+nCV2-1]); ha(curgroup+nsubgroups).XTickLabelMode='auto';
                    end
                    if curgroup == 1 
                        ylabel(ha(curgroup),ylb); ylabel(ha(curgroup+nsubgroups),ylb); 
                        xlabel(ha(curgroup+nsubgroups),sprintf('%g/%g [ %3.1f%% ] of CV_2 partitions processed',ll,nCV2,ll*100/nCV2)); 
                        legend(ha(curgroup), lg{curgroup},'Location','Best'); 
                        legend(ha(curgroup+nsubgroups), lh{curgroup},'Location','Best'); 
                    end
                    box(ha(curgroup),'on'); ha(curgroup).YGrid='on';
                    box(ha(curgroup+nsubgroups),'on'); ha(curgroup+nsubgroups).YGrid='on';
                    hold(ha(curgroup),'off');hold (ha(curgroup+nsubgroups),'off');
                end
                drawnow
            end
        end
        [RootPath, FileNames{f,d}] = fileparts(oOOCVpath);  
        ll=ll+1;
    end
end

for curclass = 1: nclass

    if inp.targscale 
        IN.minY = inp.minLbCV; IN.maxY = inp.maxLbCV; IN.revertflag = 1;
        binOOCVD{curclass} = nk_PerfScaleObj(binOOCVD{curclass}, IN);
        if LabelMode
            labelOOCV = nk_PerfScaleObj(labelOOCV, IN);
        end
    end
    if ~isempty(inp.PolyFact)
        binOOCVD{curclass} = binOOCVD{curclass} .^ (1/inp.PolyFact);
        labelOOCV = labelOOCV .^ (1/inp.PolyFact);
    end

    Results.PerformanceMeasure = char(EVALFUNC);

    switch MODEFL
        case 'classification'
            Results.BinCV2Predictions_DecisionValues{curclass}  = binOOCVD{curclass};
            Results.BinCV2Predictions_Targets{curclass}         = sign(binOOCVD{curclass});
            if size(Results.BinCV2Predictions_DecisionValues{curclass},2) > 1
                Results.MeanCV2PredictedValues{curclass}        = nm_nanmedian(Results.BinCV2Predictions_DecisionValues{curclass},2);
                Results.StdCV2PredictedValues{curclass}         = nm_nanstd(Results.BinCV2Predictions_DecisionValues{curclass},2);
            else
                Results.MeanCV2PredictedValues{curclass}        = Results.BinCV2Predictions_DecisionValues{curclass};
                Results.StdCV2PredictedValues{curclass}         = zeros(size(Results.BinCV2Predictions_DecisionValues{curclass}));
            end
            Results.CICV2PredictedValues{curclass}              = cell2mat(arrayfun( @(i) percentile(Results.BinCV2Predictions_DecisionValues{curclass}(i,:),[2.5 97.5]), ...
                                                                        1:length(Results.MeanCV2PredictedValues{curclass}),'UniformOutput',false)');
            [Results.BinProbPredictions{curclass},...
                Results.BinMajVoteProbabilities{curclass}]      = nk_MajVote(Results.BinCV2Predictions_DecisionValues{curclass},[1 -1]);
            if ~exist("indnan","var") || isempty(indnan)
                indnan = isnan(Results.MeanCV2PredictedValues{curclass});
            end
            Results.MeanCV2PredictedValues{curclass}(indnan)=nan;
            Results.StdCV2PredictedValues{curclass}(indnan)=nan;
            Results.CICV2PredictedValues{curclass}(indnan)=nan;
            Results.BinProbPredictions{curclass}(indnan)=nan;
            Results.BinMajVoteProbabilities{curclass}(indnan)=nan;

            if inp.PERM.flag==1
                binOOCVDhx = binOOCVD{curclass}(indDicho{curclass},:);
                hdx = nm_nanmedian(binOOCVDhx,2); hdx(indnan) = nan;
                PermPredPerf = zeros(1,inp.PERM.nperms);
                for curperm = 1:inp.PERM.nperms
                    PermPred = nm_nanmedian(binOOCVD_perm{curclass}(indDicho{curclass},:,curperm),2);
                    PermPredPerf(curperm) = EVALFUNC(labelDicho{curclass}(indDicho{curclass}), PermPred);
                end
                Results.PermAnal.ModelObsPerf(curclass) = EVALFUNC(labelDicho{curclass}(indDicho{curclass}), hdx);
                Results.PermAnal.ModelPermPerf{curclass} = PermPredPerf;
                Results.PermAnal.ModelPermPerfCrit{curclass} = compfun(PermPredPerf, Results.PermAnal.ModelObsPerf(curclass));
                Results.PermAnal.ModelPermSignificance(curclass) = sum(Results.PermAnal.ModelPermPerfCrit{curclass})/inp.PERM.nperms;    
            end

            if isfield(inp,'groupind')
                for g = 1:nsubgroups
                    indr = []; 
                    if iscell(inp.groupvec)
                        indg = find(strcmp(inp.groupind,inp.groupvec{g}));
                    else
                        indg = find(inp.groupind == inp.groupvec(g));
                    end
                    refstr = sprintf('%s', inp.groupnames{g});
                    if isfield(inp,'refgroup')
                        if iscell(inp.groupvec)
                            indr = find(strcmp(inp.groupind,inp.refgroup));
                        else
                            indr = find(inp.groupind == inp.refgroup);
                        end
                        refstr = sprintf('%s vs %s', inp.groupnames{g}, inp.groupnames{inp.refgroup});
                    end
                    idx = [indg; indr];
                    if isfield(inp,'groupnames'), Results.Group{g}.GroupName = refstr; end
                    Results.Group{g}.ObservedValues{curclass} = labelDicho{curclass}(idx);
                    Results.Group{g}.MeanCV2PredictedValues{curclass} = Results.MeanCV2PredictedValues{curclass}(idx);
                    Results.Group{g}.StdCV2PredictedValues{curclass}  = Results.StdCV2PredictedValues{curclass}(idx);
                    Results.Group{g}.CICV2PredictedValues{curclass}   = cell2mat(arrayfun( @(i) percentile(Results.BinCV2Predictions_DecisionValues{curclass}(idx(i),:),[2.5 97.5]), ...
                            1:numel(idx),'UniformOutput',false)');
                    if LabelMode
                        Results.Group{g}.PredictionPerformance = ALLPARAM(Results.Group{g}.ObservedValues{curclass}, Results.Group{g}.MeanCV2PredictedValues{curclass}); 
                    end
                    if inp.PERM.flag==1
                        hdx = Results.Group{g}.MeanCV2PredictedValues{curclass}; 
                        PermPredPerf = zeros(1,inp.PERM.nperms);
                        for curperm = 1:inp.PERM.nperms
                            PermPred = nm_nanmedian(binOOCVD_perm{curclass}(idx,:,curperm),2);
                            PermPredPerf(curperm) = EVALFUNC(Results.Group{g}.ObservedValues{curclass}, PermPred);
                        end
                        Results.Group{g}.PermAnal.ModelObsPerf(curclass) = EVALFUNC(Results.Group{g}.ObservedValues{curclass}, hdx);
                        Results.Group{g}.PermAnal.ModelPermPerf{curclass} = PermPredPerf;
                        Results.Group{g}.PermAnal.ModelPermPerfCrit{curclass} = compfun(PermPredPerf, Results.Group{g}.PermAnal.ModelObsPerf(curclass));
                        Results.Group{g}.PermAnal.ModelPermSignificance(curclass) = sum(Results.Group{g}.PermAnal.ModelPermPerfCrit{curclass})/inp.PERM.nperms;  
                    end
                end
            end

        case 'regression'
            
            Results.CV2PredictedValues      = binOOCVD{1};
            Results.MeanCV2PredictedValues  = nm_nanmedian(Results.CV2PredictedValues,2);
            Results.StdCV2PredictedValues   = nm_nanstd(Results.CV2PredictedValues,2);
            Results.CICV2PredictedValues    = cell2mat(arrayfun( @(i) percentile(Results.CV2PredictedValues(i,:),[2.5 97.5]), ...
                1:size(Results.CV2PredictedValues,1),'UniformOutput',false)');
            Results.MeanCV2PredictedValues(indnan)=nan;
            Results.StdCV2PredictedValues(indnan)=nan;
            Results.CICV2PredictedValues(indnan)=nan;

            Results.ErrCV2PredictedValues = Results.MeanCV2PredictedValues - labelOOCV;
            if LabelMode
                Results.Regr = nk_ComputeEnsembleProbability(Results.MeanCV2PredictedValues , labelOOCV, 1);
            else
                Results.Regr = nk_ComputeEnsembleProbability(Results.MeanCV2PredictedValues , [], 1);
            end
            if inp.PERM.flag==1
                binOOCVDhx = binOOCVD{1};
                hdx = nm_nanmedian(binOOCVDhx,2); hdx(indnan) = nan;
                PermPredPerf = zeros(1,inp.PERM.nperms);
                for curperm = 1:inp.PERM.nperms
                    PermPred = nm_nanmedian(binOOCVD_perm{1}(:,:,curperm),2);
                    PermPredPerf(curperm) = EVALFUNC(labelOOCV, PermPred);
                end
                Results.PermAnal.ModelObsPerf = EVALFUNC(labelOOCV, hdx);
                Results.PermAnal.ModelPermPerf = PermPredPerf;
                Results.PermAnal.ModelPermPerfCrit = compfun(PermPredPerf, Results.PermAnal.ModelObsPerf);
                Results.PermAnal.ModelPermSignificance = sum(Results.PermAnal.ModelPermPerfCrit)/inp.PERM.nperms;    
            end

            if nsubgroups > 1 && isfield(inp,'groupind')
                try
                    [Results.GroupComp.P, Results.GroupComp.AnovaTab, Results.GroupComp.Stats] = ...
                        anova1(Results.ErrCV2PedictedValues,inp.groupind);
                    Results.GroupComp.MultCompare = multcompare(Results.GroupComp.Stats);
                catch
                    warning('Group comparison statistics not supported!')
                end
                vec = unique(inp.groupind);
                for g = 1:numel(vec)
                    if iscell(vec)
                        indg = find(strcmp(inp.groupind, vec{g}));
                    else
                        indg = find(inp.groupind == vec(g));
                    end
                    Results.Group{g}.ObservedValues = labelOOCV(indg);
                    Results.Group{g}.Index = indg;
                    Results.Group{g}.MeanCV2PredictedValues = Results.MeanCV2PredictedValues(indg);
                    Results.Group{g}.StdCV2PredictedValues = Results.StdCV2PredictedValues(indg);
                    Results.Group{g}.CICV2PredictedValues   = cell2mat(arrayfun( @(i) percentile(Results.CV2PredictedValues(indg(i),:),[2.5 97.5]), 1:numel(indg),'UniformOutput',false)');
                    Results.Group{g}.ErrCV2PedictedValues = Results.ErrCV2PredictedValues(indg);
                    if isfield(inp,'groupnames'), Results.Group{g}.GroupName = inp.groupnames{g}; end
                    if numel(Results.Group{g}.ObservedValues)>2
                        Results.Group{g}.PredictionPerformance = EVALFUNC(Results.Group{g}.ObservedValues, Results.Group{g}.MeanCV2PredictedValues);
                        Results.Group{g}.CorrPredictObserved = corrcoef(Results.Group{g}.ObservedValues,Results.Group{g}.MeanCV2PredictedValues);
                        Results.Group{g}.CorrPredictObserved = Results.Group{g}.CorrPredictObserved(2);
                        Results.Group{g}.Regr = nk_ComputeEnsembleProbability(Results.Group{g}.MeanCV2PredictedValues, Results.Group{g}.ObservedValues, 1);
                        if inp.PERM.flag==1
                            hdx = Results.Group{g}.MeanCV2PredictedValues; 
                            PermPredPerf = zeros(1,inp.PERM.nperms);
                            for curperm = 1:inp.PERM.nperms
                                PermPred = nm_nanmedian(binOOCVD_perm{1}(indg,:,curperm),2);
                                PermPredPerf(curperm) = EVALFUNC(Results.Group{g}.ObservedValues, PermPred);
                            end
                            Results.Group{g}.PermAnal.ModelObsPerf = EVALFUNC(Results.Group{g}.ObservedValues, hdx);
                            Results.Group{g}.PermAnal.ModelPermPerf = PermPredPerf;
                            Results.Group{g}.PermAnal.ModelPermPerfCrit = compfun(PermPredPerf, Results.Group{g}.PermAnal.ModelObsPerf);
                            Results.Group{g}.PermAnal.ModelPermSignificance = sum(Results.Group{g}.PermAnal.ModelPermPerfCrit)/inp.PERM.nperms;  
                        end
                    else
                        fprintf('Need more than 2 subjects to compute meaningful performance metric.');
                    end
                end
            end
    end
end

if MULTI.flag && multiflag 
    labels = []; if LabelMode, labels = labelOOCV; end
    % Compute multiclass prediction results for current OOCV data
    Results = nk_MultiPerfComp(Results, Results.MultiCV2PredictionsLL, labels, ngroups, 'pred');
    % Permutation enabled?
    if inp.PERM.flag == 1
        PermPredPerf = zeros(1,inp.PERM.nperms);
        for curperm = 1:inp.PERM.nperms
            nDicho = 0; nDichoH = zeros(nclass+1,1); nDichoH(1) = 1;
            for curclass = 1 : nclass
                nDicho = nDicho + size(binOOCVD{curclass}(:,:,curperm),2);
                nDichoH(curclass+1) = 1 + nDicho;
            end
            multiOOCV_perm  = zeros(inp.nOOCVsubj, nDicho);
            multiClass_perm = zeros(1, nDicho);
            for curclass = 1 : nclass
                multiOOCV_perm( : , nDichoH(curclass) : nDichoH(curclass+1) - 1) = binOOCVD{curclass}(:,:,curperm);
                multiClass_perm(1, nDichoH(curclass) : nDichoH(curclass+1) - 1) = ones(1,size(binOOCVD{curclass}(:,:,curperm),2)) * curclass;
            end
            PermPredPerf(curperm) = nk_MultiEnsPerf(multiOOCV_perm, sign(multiOOCV_perm), labels, multiClass_perm, ngroups);
        end
        Results.MultiClass.PermAnal.ModelObsPerf = Results.MultiClass.performance.BAC_Mean;
        Results.MultiClass.PermAnal.ModelPermPerf = PermPredPerf;
        Results.MultiClass.PermAnal.ModelPermPerfCrit = compfun(PermPredPerf, Results.MultiClass.PermAnal.ModelObsPerf);
        Results.MultiClass.PermAnal.ModelPermSignificance = sum(Results.MultiClass.PermAnal.ModelPermPerfCrit)/inp.PERM.nperms;  
    end
end

