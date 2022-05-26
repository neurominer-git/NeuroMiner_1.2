function [Results, FileNames, RootPath] = nk_MLInterpreter(inp)
% =========================================================================
% FORMAT Results = nk_MLInterpreter(inp)
% =========================================================================
% Interpretable ML module
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, last modified 08/2020

global SVM RFE MULTI MODEFL CV EVALFUNC SCALE SAV CVPOS 
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
ngroups         = inp.ngroups;
analysis        = inp.analysis;
GridAct         = inp.GridAct;
batchflag       = inp.batchflag;
algostr         = GetMLType(SVM);
[ylm, ylb]      = nk_GetScaleYAxisLabel(SVM);

% Setup CV2 and CV1 counters and data containers:
[ix, jx]        = size(CV.TrainInd);
[iy, jy]        = size(CV.cvin{1,1}.TrainInd);
interpAll       = cell(nclass,1);
ll              = 1; 
totLearn        = 0;

if ~exist('GridAct','var') || isempty(GridAct), GridAct = nk_CVGridSelector(ix,jx); end
if ~exist('batchflag','var') || isempty(batchflag), batchflag = false; end

CL = getNMcolors;

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
        predCV2 = cell(nclass,1);
        switch inp.method
            case 'posneg'
                interpCV2 = cell(nclass,2);
            case 'median'
                interpCV2 = cell(nclass,1);
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
                            
                            %% Step 2: Apply trained model to OOCV data 
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
                                    [ TR , CV1, CV2 ] = nk_ReturnAtOptPos(mapY.Tr{k,l}{hix},  mapY.CV{k,l}{hix}, mapY.Ts{k,l}{hix}, mapYocv.Ts{k,l}{hix}, Param{1}(k,l,hix), pnt);                                       
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

                                    uD = zeros(size(TsInd,1),ul);

                                    % Loop through feature subspaces
                                    for u=1:ul
                                        % Extract features according to mask
                                        TR_star = nk_ExtractFeatures(TR, F, [], u);
                                        CV2_star = nk_ExtractFeatures(CV2, F, [], u);
                                        if ~fndMD
                                            fprintf('Compute optimal model...');
                                            [~, MD{h}{m}{k,l}{u}] = nk_GetParam2(TR_star, modelTrL, sPs, 1);
                                            fprintf(' done.');
                                        end
                                        fprintf(' Predict CV2 test data...')
                                        % Apply model to CV2 test data 
                                        [~, ~, uD(:,u)] = nk_GetTestPerf(CV2_star, ones(size(CV2_star,1),1), F(:,u), MD{h}{m}{k,l}{u}, Ymodel, 1);
                                        fprintf(' done.')
                                    end
                                    
                                    if detrendfl 
                                        fprintf(' Detrend CV2 predictions.')
                                        beta = GD.Detrend{Pspos(m)}.beta;
                                        p = GD.Detrend{Pspos(m)}.p;
                                        uD(:,u) = nk_DetrendPredictions2(beta, p, uD(:,u)); 
                                    end    
                                    predCV2{h} = [predCV2{h} uD];
                                end
                            end
                        end
                    end

                    for h=1:nclass  % Loop through binary comparisons

                        if nclass > 1, fprintf('\n\n*** %s #%g ***',algostr, h); end

                        %% Step 1: Get optimal model parameters
                        % Retrieve optimal parameters from precomputed analysis structure
                        % Differentiate according to binary or multi-group mode
                        [~, ~, nP, ~] = nk_GetModelParams2(analysis, multiflag, ll, h, inp.curlabel);
                        
                        if ~fndMD , MD{h} = cell(nP,1); end
                        predOCV = cell(1,nclass);

                        for q=1:numel(tInd) % Loop through CV2/OOCV cases 
                            
                            for nx = 1:numel(inp.X)
                                % Create artificial data and add it as "OOCV
                                % data" (inp.X(nx).Yocv) to the input structure. 
                                % if OOCV data has been selected for interpretation, then create
                                % inp.X(nx).Yocv2 and deal with this scenario in nk_ApplyTrainedPreproc 
                                % (see there)
                                inp = nk_CreateData4MLInterpreter(inp, nx, TrInd, tInd(q,:));
                            end
                            switch inp.method
                                case 'posneg'
                                    inp.desc_oocv{1} = sprintf('%g%%-percentile modification', inp.upper_thresh);
                                    inp.desc_oocv{2} = sprintf('%g%%-percentile modification', inp.lower_thresh);
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

                                        % get training data using pointers
                                        % Either (a) only CV1 training data, (b) CV1
                                        % training and test data, (c) CV1 training & test
                                        % data as well as CV2 test data              
                                        if BINMOD, hix = h; else, hix = 1; end
                                        [ TR , ~, ~, OCV ] = nk_ReturnAtOptPos(mapY.Tr{k,l}{hix},  mapY.CV{k,l}{hix}, mapY.Ts{k,l}{hix}, mapYocv.Ts{k,l}{hix}, Param{1}(k,l,hix), pnt);    
                                        switch inp.method
                                            case 'posneg'
                                                % Loop through feature subspaces
                                                for u=1:ul
            
                                                    % Extract features according to mask
                                                    TR_star   = nk_ExtractFeatures(TR, F, [], u);
                                                    OCV_star_pos  = nk_ExtractFeatures(OCV{1}, F, [], u);
                                                    OCV_star_neg  = nk_ExtractFeatures(OCV{2}, F, [], u);
                                                    % Apply trained model to
                                                    % artificial data and generate
                                                    % predictions for later
                                                    % evaluation
                                                    [~, ~, uD_pos(:,u)] = nk_GetTestPerf(OCV_star_pos, ones(size(OCV_star_pos,1),1), F(:,u), MD{h}{m}{k,l}{u}, TR_star, 1);
                                                    [~, ~, uD_neg(:,u)] = nk_GetTestPerf(OCV_star_neg, ones(size(OCV_star_neg,1),1), F(:,u), MD{h}{m}{k,l}{u}, TR_star, 1);
        
                                                    % Detrend regressor predictions, if
                                                    % required by user input
                                                    if detrendfl 
                                                        beta = GD.Detrend{Pspos(m)}.beta;
                                                        p = GD.Detrend{Pspos(m)}.p;
                                                        uD_pos(:,u) = nk_DetrendPredictions2(beta, p, uD_pos(:,u));
                                                        uD_neg(:,u) = nk_DetrendPredictions2(beta, p, uD_neg(:,u));
                                                    end    
                                                end
                                                predOCV_pos{h} = [predOCV_pos{h} uD_pos];
                                                predOCV_neg{h} = [predOCV_neg{h} uD_neg];

                                            case 'median'

                                                % Loop through feature subspaces
                                                for u=1:ul
            
                                                    % Extract features according to mask
                                                    TR_star   = nk_ExtractFeatures(TR, F, [], u);
                                                    OCV_star  = nk_ExtractFeatures(OCV, F, [], u);
                                                    % Apply trained model to
                                                    % artificial data and generate
                                                    % predictions for later
                                                    % evaluation
                                                    [~, ~, uD(:,u)] = nk_GetTestPerf(OCV_star, ones(size(OCV_star,1),1), F(:,u), MD{h}{m}{k,l}{u}, TR_star, 1);
        
                                                    % Detrend regressor predictions, if
                                                    % required by user input
                                                    if detrendfl 
                                                        beta = GD.Detrend{Pspos(m)}.beta;
                                                        p = GD.Detrend{Pspos(m)}.p;
                                                        uD(:,u) = nk_DetrendPredictions2(beta, p, uD(:,u)); 
                                                    end    
                                                end
                                                predOCV{h} = [predOCV{h} uD];
                                        end
                                        
                                    end
                                end
                            end
                        end

                    end

                    fprintf('\nSaving %s', oMLIpath); save(oMLIpath,'binOOCVDh','operm','ofold');
                    if saveparam, fprintf('\nSaving %s', OptModelPath); save(OptModelPath, 'MD', 'ofold','operm'); end
                end
                
            case 1
                
                vpth = deblank(matfiles{f,d});
                if isempty(vpth) || ~exist(vpth,'file') && GridAct(f,d)
                    error(['No valid OOCVdatamat detected for CV2 partition ' '[' num2str(f) ', ' num2str(d) ']!']);
                else
                    [~,vnam] = fileparts(vpth);
                    fprintf('\n\nLoading ML interpretation results for CV2 partition [ %g, %g ]:', f, d);
                    fprintf('\n%s',vnam);
                    load(vpth)
                end 
        end
        
        for h=1:nclass
            if meanflag == 2
                fprintf('\nCompute mean of base learners'' outputs of current CV2 partition and add them to the ensemble matrix.')
                binOOCVD{h} = [binOOCVD{h} nm_nanmedian(binOOCVDh{h},2)];
            else
                fprintf('\nAdd all base learners'' outputs to the ensemble matrix without averaging.')
                binOOCVD{h} = [binOOCVD{h} binOOCVDh{h}];
            end
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

            if LabelMode,    
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
                    binOOCVDhx_ll = binOOCVDh{curclass}(indDicho{curclass},:);
                    binOOCVDhx = binOOCVD{curclass}(indDicho{curclass},:);
                    hrx_ll = sign(nm_nansum(sign(binOOCVDhx_ll),2)); 
                    hdx_ll = nm_nanmedian(binOOCVDhx_ll,2); hdx_ll(indnan) = nan;
                    hrx = sign(nm_nansum(sign(binOOCVDhx),2)); 
                    hdx = nm_nanmedian(binOOCVDhx,2); hdx(indnan) = nan;
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
                end
            else
                hdx_ll = nm_nanmedian(binOOCVDh{1},2); hdx_ll(indnan) = nan;
                hdx = nm_nanmedian(binOOCVD{1},2); hdx(indnan) = nan;
                if ll==1, Results.RegrLabels = labelOOCV; Results.RegrLabels(indnan) = nan; end
                Results.CV2Performance_PredictedValues(ll) = EVALFUNC(Results.RegrLabels, hdx_ll);
                Results.CV2Performance_PredictedValues_History(ll) = EVALFUNC(Results.RegrLabels, hdx);
            end
            % Disply progress information in figure
            if ~batchflag
                if ~exist('hu','var') || isempty(hu) 
                    hu = findobj('Tag','OOCV');
                    if isempty(hu), 
                        sz = get(0,'ScreenSize');
                        win_wdth = sz(3)/3; win_hght = sz(4)/1.25; win_x = sz(3)/2 - win_wdth/2; win_y = sz(4)/2 - win_hght/2;
                        hu = figure('Tag','OOCV', ...
                            'NumberTitle','off', ...
                             'MenuBar','none', ...
                             'Color', [0.9 0.9 0.9], ...
                             'Position', [win_x win_y win_wdth win_hght]); 
                    end
                end
                hu.Name = sprintf('OOCV prediction viewer: %s',inp.analysis_id);
                set(0,'CurrentFigure',hu); clf; ha = tight_subplot(2,1,0.03,[.1 .05],[.1 .05]); hold(ha(1),'on');hold(ha(2),'on');
                lg = []; hc = []; lh = []; hd = [];
                for curclass=1:nclass
                    switch MODEFL
                        case 'classification'
                            hc(end+1) = plot(ha(1),Results.BinCV2Performance_Targets_History(curclass,:),CLtarg_m,'Color', CL(curclass,:));
                            hc(end+1) = plot(ha(1),Results.BinCV2Performance_DecisionValues_History(curclass,:), CLdec_m,'Color', CL(curclass,:));
                            lg{end+1} = sprintf('Classifier %g: prediction target performance of overall ensemble.',curclass); 
                            lg{end+1} = sprintf('Classifier %g: prediction score performance of overall ensemble.',curclass); 
                            hd(end+1) = plot(ha(2),Results.BinCV2Performance_Targets(curclass,:),CLtarg_m,'Color', CL(curclass,:));
                            hd(end+1) = plot(ha(2),Results.BinCV2Performance_DecisionValues(curclass,:), CLdec_m,'Color', CL(curclass,:));
                            lh{end+1} = sprintf('Classifier %g: prediction target performance of current CV2 ensemble.',curclass); 
                            lh{end+1} = sprintf('Classifier %g: prediction score performance of current CV2 ensemble.',curclass); 
                        case 'regression'
                            hd(end+1) = plot(ha(1),Results.CV2Performance_PredictedValues_History, CLtarg_m, 'Color', CL(curclass,:));
                            lg{end+1} = sprintf('Predictor performance of current ensemble.'); 
                            hd(end+1) = plot(ha(2),Results.CV2Performance_PredictedValues, CLtarg_m, 'Color', CL(curclass,:));
                            lh{end+1} = sprintf('Predictor performance of current ensemble.'); 
                    end

                end
                if MULTI.flag && multiflag,plot(Results.MultiCV2Performance_History,'k+'); end
                ylim(ha(1),ylm); ylabel(ha(1),ylb); ha(1).YTickLabelMode='auto';
                ylim(ha(2),ylm); ylabel(ha(2),ylb); ha(2).YTickLabelMode='auto';
                if ll_start+nCV2-1 > ll_start
                    xlim(ha(1),[ll_start ll_start+nCV2-1]); ha(1).XTickLabelMode='auto';
                    xlim(ha(2),[ll_start ll_start+nCV2-1]); ha(2).XTickLabelMode='auto';
                end
                
                xlabel(ha(2),sprintf('%g/%g [ %3.1f%% ] of CV_2 partitions processed',ll,nCV2,ll*100/nCV2)); 
                legend(ha(1), lg,'Location','Best'); 
                legend(ha(2), lh,'Location','Best'); 
                box(ha(1),'on'); ha(1).YGrid='on';
                box(ha(2),'on'); ha(2).YGrid='on';
                drawnow      
            end
        end
        [RootPath, FileNames{f,d}] = fileparts(oMLIpath);  
        ll=ll+1;
    end
end

for curclass = 1: nclass

    if inp.targscale
        IN.minY = inp.minLbCV; IN.maxY = inp.maxLbCV; IN.revertflag = 1;
        binOOCVD{curclass} = nk_PerfScaleObj(binOOCVD{curclass}, IN);
        labelOOCV = nk_PerfScaleObj(labelOOCV, IN);
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
            Results.MeanCV2PredictedValues{curclass}            = nm_nanmedian(Results.BinCV2Predictions_DecisionValues{curclass},2);
            Results.StdCV2PredictedValues{curclass}             = nm_nanstd(Results.BinCV2Predictions_DecisionValues{curclass},2);
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

            if isfield(inp,'groupind')
                vec = unique(inp.groupind);
                for g = 1:numel(vec)
                    indg = find(inp.groupind == vec(g));
                    if LabelMode
                        if curclass == 1
                            Results.Group{g}.ObservedValues{curclass} = labelDicho{curclass}(indg);
                        end
                    end
                    if isfield(inp,'groupnames'), Results.Group{g}.GroupName = inp.groupnames(g); end
                    Results.Group{g}.MeanCV2PredictedValues{curclass}   = Results.MeanCV2PredictedValues{curclass}(indg);
                    Results.Group{g}.StdCV2PredictedValues{curclass}    = Results.StdCV2PredictedValues{curclass}(indg);
                    Results.Group{g}.CICV2PredictedValues{curclass}     = cell2mat(arrayfun( @(i) percentile(Results.BinCV2Predictions_DecisionValues{curclass}(indg(i),:),[2.5 97.5]), ...
                            1:numel(indg),'UniformOutput',false)');
                    if LabelMode, Results.Group{g}.PredictionPerformance= ALLPARAM(Results.Group{g}.ObservedValues{curclass}, Results.Group{g}.MeanCV2PredictedValues{curclass}); end
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
            if inp.ngroups > 1 && isfield(inp,'groupind')
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
                    Results.Group{g}.MeanCV2PredictedValues = Results.MeanCV2PredictedValues(indg);
                    Results.Group{g}.StdCV2PredictedValues = Results.StdCV2PredictedValues(indg);
                    Results.Group{g}.CICV2PredictedValues   = cell2mat(arrayfun( @(i) percentile(Results.CV2PredictedValues(indg(i),:),[2.5 97.5]), 1:numel(indg),'UniformOutput',false)');
                    Results.Group{g}.ErrCV2PedictedValues = Results.ErrCV2PredictedValues(indg);
                    if isfield(inp,'groupnames'), Results.Group{g}.GroupName = inp.groupnames(g); end
                    if numel(Results.Group{g}.ObservedValues)>2
                        Results.Group{g}.PredictionPerformance = EVALFUNC(Results.Group{g}.ObservedValues, Results.Group{g}.MeanCV2PredictedValues);
                        Results.Group{g}.CorrPredictObserved = corrcoef(Results.Group{g}.ObservedValues,Results.Group{g}.MeanCV2PredictedValues);
                        Results.Group{g}.CorrPredictObserved = Results.Group{g}.CorrPredictObserved(2);
                    else
                        fprintf('Need more than 2 subjects to compute meaningful performance metric.');
                    end
                end
            end
    end
end

if MULTI.flag && multiflag 
    labels = []; if LabelMode, labels = labelOOCV; end
    Results = nk_MultiPerfComp(Results, Results.MultiCV2PredictionsLL, labels, ngroups, 'pred');
end

