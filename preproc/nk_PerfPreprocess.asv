function [tY, Pnt, paramfl, tYocv] = nk_PerfPreprocess(Y, inp, labels, ...
                                                       paramfl, Yocv, Cocv)
% =========================================================================
% [tY, Pnt, paramfl, tYocv] = nk_PerfPreprocess(Y, inp, labels, ...
%                                                      paramfl, Yocv, Cocv)
% =========================================================================
% Core function of NM's preprocessing module that can be run in training
% mode as well as in OOCV mode and executes the sequence of preprocessing
% steps on the data as defined by the user's configuration.
%
% INPUTS:
% -------
% Y         = data matrix as [m x n] double, m = samples, n = dimensions
% inp 		= Parameter structure for the computational part
% labels    = n x m label vector matrix
% paramfl   = Parameter structure for the script execution part
% Yocv      = Independent test data
% Cocv      = Calibration data [currently not used]
% 
% OUTPUTS:
% --------
% tY        = the preprocessed data
% Pnt       = trained preprocessing parameters/models
% paramfl   = modified script execution parameters
% tYocv     = the preprocessed independent test data
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 07/2022

global PREPROC MODEFL MULTI CV xCV RAND VERBOSE TEMPL SVM MULTILABEL CALIB

% Initialize runtime variables
i       = inp.f; % Curration permutation
j       = inp.d; % Current fold
kbin    = inp.nclass;

% Check whether the function is used by the simulation module or not and
% transfer the right CV structure to tCV
if isfield(inp,'simFlag') && inp.simFlag
    tCV = xCV;
else
    tCV = CV;
end

[iy,jy] = size(tCV(1).cvin{i,j}.TrainInd); 
TrInd   = tCV(1).TrainInd{i,j}; 
TsInd   = tCV(1).TestInd{i,j}; 
TsI     = cell(kbin,1);

if ~exist('paramfl','var'), paramfl.use_exist = false; end
cv2flag = false; if isfield(PREPROC,'CV2flag') && (PREPROC.CV2flag - 1) == 1; cv2flag = true; end
tYocv   = []; if exist('Yocv','var') && ~isempty(Yocv), tYocv.Ts = cell(iy,jy); else, Yocv = []; end
if ~exist('Cocv','var'), Cocv = []; end

if isfield(inp, 'C')
    Cocv = inp.C{1,1}.Y;
    CALIB = inp.C{1,1};
end

% Set binary/regression or multi-group processing mode
if iscell(PREPROC)
    BINMOD = PREPROC{1}.BINMOD;
else
    BINMOD = PREPROC.BINMOD;
end
if BINMOD || strcmp(MODEFL,'regression')
    ukbin = kbin;   if VERBOSE; fprintf('\nProcessing Mode: binary / regression preprocessing'); end
else
    ukbin = 1;      
    if VERBOSE; fprintf('\nProcessing Mode: multi-group preprocessing'); end
end

if VERBOSE
    if iscell(Y)
        fprintf('\nMultiple shelfs of input data detected')
        for ii=1:size(Y,1)
            fprintf('\nShelf [ %2g ]: Original dimensionality: %g', ii, size(Y{ii},2)); 
        end
    else
        fprintf('\nOriginal dimensionality: %g', size(Y,2)); 
    end    
end

% Set up data containers
tY.Tr = cell(iy,jy);
tY.CV = cell(iy,jy);
tY.Ts = cell(iy,jy);
    
% Labels & Indices
tY.TrL = cell(iy,jy);
tY.CVL = cell(iy,jy);
tY.TrInd = cell(iy,jy);
tY.CVInd = cell(iy,jy);
multoocv = false; if iscell(Yocv) && numel(Yocv)>1, multoocv = true; end

% Eventually, apply spatial operations to image data
% (This function will be extended beyond smoothing ops on nifti data)
[sY, sYocv, sCocv, inp, optfl, ukbin, uBINMOD, BINMOD] = ...
    nk_PerfPreprocessSpatial( Y, Yocv, Cocv, inp, paramfl, BINMOD, kbin, ukbin);

if ~BINMOD && isfield(paramfl,'PXopt') && numel(paramfl.PXopt)>1
    % Here, we force a multi-group processing mode but map the multi-group processed data into binary containers
    uBINMOD = 0; 
    ukbin = 1; 
    if VERBOSE; fprintf('\nProcessing Mode: multi-group preprocessing, but no multi-group classifier requested'); end
end

% Generate template parameters (e.g. for Procrustes rotation)
if isfield(paramfl,'templateflag') && paramfl.templateflag 
    TEMPL = nk_GenTemplParam(PREPROC, tCV, MODEFL, RAND, sY, inp, kbin);
else
    clear TEMPL
end

% The order of statements here is critical: TsI should receive the original
% TsInd before it is modified below!
if ~isempty(MULTI) && MULTI.flag, tY.mTsL = labels(TsInd,:); end
for u=1:ukbin, TsI{u} = TsInd; end

if MULTILABEL.flag && size(labels,2)>1
    lb = MULTILABEL.sel;
else
    lb = 1;
end

for u=1:kbin
    switch MODEFL
        case 'regression'
            TsL     = labels(TsInd,lb);
            TsInd   = true(size(TsL,1),1);
        case 'classification'
            if RAND.Decompose == 9
                TsL     = labels(TsInd,lb);
                TsInd   = true(size(TsL,1),1);
            else
                TsL     = tCV.classnew{i,j}{u}.label(:,lb);
                TsInd   = tCV.classnew{i,j}{u}.ind;
            end
    end
    tY.TsL{u} = TsL;
    tY.TsInd{u} = TsInd;
end

if ~exist('paramfl','var') || ~paramfl.found || ~paramfl.use_exist
    Pnts = struct('data_ind', [], ...
                  'train_ind', [], ...
                  'nP', [], ...
                  'nA', [], ...
                  'TrainedParam', []);
    Pnt = repmat(Pnts,iy,jy,ukbin);
else
    Pnt = paramfl.Param;
    paramfl.Param = [];
end

% train parameters eventually only for a subset of partitions
if isfield(inp,'CV1p')
    sta_iy = inp.CV1p(1); stp_iy = inp.CV1p(2);
    sta_jy = inp.CV1f(1); stp_jy = inp.CV1f(2);
else
    sta_iy = 1; stp_iy = iy;
    sta_jy = 1; stp_jy = jy;
end
oocvonly = false;
if isfield(inp,'OOCVonly') && inp.OOCVonly, oocvonly= true; end

for k=sta_iy:stp_iy % Inner permutation loop

    for l=sta_jy:stp_jy % Inner CV fold loop
        
        tElapsed = tic;
        fprintf('\nWorking on CV1 [%2g, %2g ]: Prepare data', k, l);
        
        for u=1:ukbin % Binary comparison loop depending on PREPROC.BINMOD & FBINMOD
            
            if ischar(Pnt(k,l,u).TrainedParam) && exist(Pnt(k,l,u).TrainedParam,'file')
                load( Pnt(k,l,u).TrainedParam );
            else
                oTrainedParam = Pnt(k,l,u).TrainedParam;
            end
            
            %Check whether optimised parameter space exists
            PREPROC = nk_SetParamChain(paramfl, u, PREPROC);

            if optfl
                usY = sY{u}; 
                if ~isempty(Yocv)
                    if multoocv
                        usYocv = cell(1,numel(sYocv));
                        for n=1:numel(sYocv)
                            usYocv{n} = sYocv{n}{u};
                        end
                    else
                        usYocv = sYocv{u}; 
                    end
                end
                if ~isempty(Cocv) 
                    usCocv = sCocv{u}; 
                end
                if isfield(inp,'Yw'), InputParam.Yw = inp.Yw{u}; end
            else
                usY = sY;
                if ~isempty(Yocv), usYocv = sYocv; end
                if ~isempty(Cocv), usCocv = sCocv; end
                if isfield(inp,'Yw'), InputParam.Yw = inp.Yw; end
            end
            paramfl.P{u} = nk_ReturnParamChain(PREPROC, 1); 
            
            %% Push data into partition / folds
            % Define pointers to data
            if ~uBINMOD || strcmp(MODEFL,'regression')
                % Multi-group label: 1, 2 ,3, ...
                TrX = TrInd(tCV.cvin{i,j}.TrainInd{k,l});
            else
                % Binary label
                TrX = TrInd(tCV.class{i,j}{u}.TrainInd{k,l});
            end
            TrLX = labels(TrX,lb);
            % Assign data to containers
            mult_contain = false; iOCV=[]; 
            if oocvonly
                TrI = TrInd(tCV.cvin{i,j}.TrainInd{k,l});
                TrL = labels(TrI,lb);
                if size(TrI,2)>2, TrI = TrI'; end
                if iscell(usY)
                    n_usY = numel(usY); mult_contain = true; 
                    vTr = cell(n_usY,1); vTs = vTr;
                    for pu = 1:n_usY
                        % Training & OOCV data
                        vTr{pu} = usY{pu}(TrX,:); 
                        if pu == 1
                            % Remove cases which are completely NaN
                            [vTr{pu}, TrLX, SrcParam.iTrX] = nk_ManageNanCases(vTr{pu}, TrLX);
                            if multoocv
                                for n=1:numel(usYocv)
                                    [vTs{pu,n}, ~, SrcParam.iOCV{n}] = nk_ManageNanCases(usYocv{n}{pu});
                                end
                            else
                                [vTs{pu}, ~, SrcParam.iOCV] = nk_ManageNanCases(usYocv{pu});
                            end
                        else
                            vTr{pu} = nk_ManageNanCases(vTr{pu});
                            if multoocv
                                for n=1:numel(usYocv)
                                    vTs{pu,n} = nk_ManageNanCases(usYocv{n}{pu});
                                end
                            else
                                vTs{pu} = nk_ManageNanCases(usYocv{pu});
                            end
                        end
                        % Calibration data
                        if ~isempty(Cocv) 
                            InputParam.C{pu} = nk_ManageNanCases(usCocv{pu}); 
                        end
                    end
                else
                    vTr = usY(TrX,:); 
                    % Remove cases which are completely NaN
                    [vTr, TrLX, SrcParam.iTrX] = nk_ManageNanCases(vTr, TrLX);
                    if multoocv
                        for n=1:numel(usYocv)
                            [vTs{n}, ~, SrcParam.iOCV{n}] = nk_ManageNanCases(usYocv{n});
                        end
                    else
                        [vTs, ~, SrcParam.iOCV] = nk_ManageNanCases(usYocv);
                    end
                    % Calibration data
                    if ~isempty(Cocv), InputParam.C = nk_ManageNanCases(usCocv); end
                end
            else
                TrI = TrInd(tCV.cvin{i,j}.TrainInd{k,l});
                CVI = TrInd(tCV.cvin{i,j}.TestInd{k,l});
                
                TrL = labels(TrI,lb);
                CVL = labels(CVI,lb);
                if size(TrI,2)>2
                    TrI = TrI';
                    CVI = CVI';
                    TsI{u} = TsI{u}';
                end
                if iscell(usY)
                    n_usY = numel(usY); mult_contain = true; 
                    vTr = cell(n_usY,1); if ~isempty(Yocv), vTs = cell(n_usY,4); else, vTs = cell(n_usY,3); end
                    for pu = 1:n_usY
                        
                        % Training & CV data
                        vTr{pu} = usY{pu}(TrX,:); vTs{pu,1} = usY{pu}(TrI,:); vTs{pu,2} = usY{pu}(CVI,:); vTs{pu,3} = usY{pu}(TsI{u},:); 
                        
                        if pu == 1
                            % Remove cases which are completely NaN
                            [vTr{pu}, TrLX, SrcParam.iTrX] = nk_ManageNanCases(vTr{pu}, TrLX);
                            [vTs{pu,1}, TrL, SrcParam.iTr] = nk_ManageNanCases(vTs{pu,1}, TrL);
                            [vTs{pu,2}, CVL, SrcParam.iCV] = nk_ManageNanCases(vTs{pu,2}, CVL);
                           
                        else
                            vTr{pu} = nk_ManageNanCases(vTr{pu});
                            vTs{pu,1} = nk_ManageNanCases(vTs{pu,1});
                            vTs{pu,2} = nk_ManageNanCases(vTs{pu,2});
                        end
                        [vTs{pu,3}, ~, SrcParam.iTs] = nk_ManageNanCases(vTs{pu,3});
                        % Independent test data ( can be multiple
                        % containers)
                        if ~isempty(Yocv)
                             if multoocv
                                 iOCV = cell(1,numel(usYocv));
                                 for n=1:numel(usYocv)
                                     [vTs{pu,3+n}, ~, iOCV{n}] = nk_ManageNanCases(usYocv{pu,2});
                                 end
                             else
                                [vTs{pu,4}, ~, iOCV] = nk_ManageNanCases(usYocv{pu});
                             end
                        end
                        % Calibration data
                        if ~isempty(Cocv) 
                            InputParam.C{pu} = nk_ManageNanCases(usCocv{pu}); 
                        end
                    end
                else
                    % Training & CV data
                    vTr = usY(TrX,:); vTs{1} = usY(TrI,:); vTs{2} = usY(CVI,:); vTs{3} = usY(TsI{u},:);
                    
                    % Remove cases which are completely NaN
                    [vTr, TrLX, SrcParam.iTrX] = nk_ManageNanCases(vTr, TrLX);
                    [vTs{1}, TrL, SrcParam.iTr] = nk_ManageNanCases(vTs{1}, TrL);
                    [vTs{2}, CVL, SrcParam.iCV] = nk_ManageNanCases(vTs{2}, CVL);
                    [vTs{3}, ~, SrcParam.iTs] = nk_ManageNanCases(vTs{3});
                    
                    % Independent test data
                    if ~isempty(Yocv)
                        if multoocv
                             iOCV = cell(1,numel(usYocv));
                             for n=1:numel(usYocv)
                                 [vTs{3+n}, ~, iOCV{n}] = nk_ManageNanCases(usYocv{n});
                             end
                        else
                            [vTs{4},~, iOCV] = nk_ManageNanCases(usYocv); 
                        end
                    end
                    % Calibration data
                    if ~isempty(Cocv), InputParam.C = nk_ManageNanCases(usCocv); end
                end
            end
            clear usYocv usCocv
            InputParam.Tr = vTr; InputParam.Ts = vTs;

            if VERBOSE; fprintf('\n-----------------------------------------------------------------------------------'); 
                switch BINMOD
                    case {1,3,4}
                        if strcmp(MODEFL,'regression') || RAND.Decompose == 9
                            fprintf('\nWorking on data partition: CV2 [%2g,%2g], CV1 [%2g,%2g]', ...
                                i, j, k, l)
                        else
                            fprintf('\nWorking on data partition: CV2 [%2g,%2g], CV1 [%2g,%2g, %s]', ...
                                i, j, k, l, CV.class{i,j}{u}.groupdesc)
                        end
                    case {0,2}
                            fprintf('\nWorking on data partition: CV2 [%2g,%2g], CV1 [%2g,%2g]', ...
                                i, j, k, l)
                end
            end
            
            %% Generate SrcParam structure
            SrcParam.oocvonly           = oocvonly;
            SrcParam.TrX                = TrX;
            if ~oocvonly
                SrcParam.TrI                = TrI;
                SrcParam.CVI                = CVI;
                SrcParam.TsI                = TsI{u};
            end
            SrcParam.u                  = u;
            SrcParam.binmult            = BINMOD;
            SrcParam.CV1perm            = k;
            SrcParam.CV1fold            = l;
            SrcParam.covars             = inp.covars;
            if isfield(inp,'covars_rep')
                SrcParam.covars_oocv        = inp.covars_rep;
            elseif isfield(inp,'covars_oocv')
                SrcParam.covars_oocv        = inp.covars_oocv;
            else 
                inp.covars_oocv = [];
            end
            if isfield(inp,'C') && isfield(inp.C{1,1}, 'covars')
                SrcParam.covars_cocv        = inp.C{1,1}.covars;
            else 
                inp.covars_cocv = [];
            end
            SrcParam.iOCV               = iOCV; % To resolve bug in nk_GenPreprocSequence.m reported by Mark Dong (29/09/2021)

            % Run ADASYN if needed
            if isfield(SVM,'ADASYN') && SVM.ADASYN.flag == 1
                if VERBOSE, fprintf('\nUsing ADASYN to generate synthetic training data for partition CV2 [%g, %g], CV [%g, %g]', i, j, k, l); else; fprintf('\t...ADASYN'); end
                Covs = [];
                % Do we have covars? if so, they have to be integrated
                % into the creation of synthetic data.
                if ~isempty(SrcParam.covars)
                    Covs = SrcParam.covars( SrcParam.TrX(~SrcParam.iTrX),:);
                    % Check if covariance matrix contains NaNs and
                    % impute them.
                    if any(isnan(Covs(:)))
                        IN = struct('method','seuclidean','k',7,'X',Covs);
                        Covs = nk_PerfImputeObj(Covs, IN);
                    end
                end
                if mult_contain
                    vTrSyn = cell(n_usY,1); LabelSyn = cell(n_usY,1); CovarsSyn = cell(n_usY,1); 
                    for pu=1:n_usY
                        [ vTrSyn{pu}, LabelSyn{pu}, CovarsSyn{pu}, SrcParam.adasynused(pu) ] = nk_PerfADASYN( vTr{pu}, TrLX, SVM.ADASYN, Covs, true );
                    end
                else
                   [ vTrSyn, LabelSyn, CovarsSyn, SrcParam.adasynused ] = nk_PerfADASYN( vTr, TrLX, SVM.ADASYN, Covs, true);
                end 
                InputParam.TrSyn = vTrSyn;
                SrcParam.TrainLabelSyn = LabelSyn;
                if ~isempty(SrcParam.covars), SrcParam.covarsSyn = CovarsSyn; end
            end
            
            switch MODEFL
                case 'classification'
                    switch RAND.Decompose 
                        case 1
                            SrcParam.BinaryTrainLabel   = [ones(sum(TrL == tCV.class{i,j}{u}.groups(1)),1); -1*ones(sum(TrL == tCV.class{i,j}{u}.groups(2)),1)];
                            SrcParam.BinaryCVLabel      = [ones(sum(CVL == tCV.class{i,j}{u}.groups(1)),1); -1*ones(sum(CVL == tCV.class{i,j}{u}.groups(2)),1)];
                        case 2
                            SrcParam.BinaryTrainLabel  = TrL ~=0; SrcParam.BinaryTrainLabel (TrL == 0) =-1;
                            SrcParam.BinaryCVLabel  = CVL ~=0; SrcParam.BinaryCVLabel (TrL == 0) =-1;
                    end
                    SrcParam.MultiTrainLabel    = TrL;
                    SrcParam.MultiCVLabel       = CVL;
                case 'regression'
                    SrcParam.TrainLabel         = TrL;
                    SrcParam.CVLabel            = CVL;
            end

            %% Generate and execute for given CV1 partition preprocessing sequence
            [InputParam, oTrainedParam, SrcParam] = nk_GenPreprocSequence(InputParam, PREPROC, SrcParam, oTrainedParam);
            
            %% Write preprocessed data to mapY structure
            if isfield(paramfl,'PREPROC') && isfield(paramfl,'PXfull') && ~isempty(paramfl.PXopt{u})
                % Here an optimized parameter space exists that has
                % been used to limit the computation to the unique
                % parameter cominations. We create and store the pointers
                % and used them later on to retrieve the right preproc
                % version of the data and the respective preproc params
                [ Pnt(k,l,u).data_ind, ...
                  Pnt(k,l,u).train_ind, ...
                  Pnt(k,l,u).nP, ...
                  Pnt(k,l,u).nA] = nk_ParamReplicator(paramfl.P{u}, paramfl.PXopt{u}, paramfl.PREPROC, numel(oTrainedParam));
            end
            if oocvonly
                if multoocv
                    ocv = [];
                    for n=1:numel(Yocv)
                        ocv = [ocv squeeze(InputParam.Ts(:,n,:))];
                    end
                else
                    ocv = squeeze(InputParam.Ts(:,1,:));
                end
                if  any(any(cellfun(@isempty,trd))) 
                    error('Empty data matrices returned by preprocessing pipeline. Check your settings')
                end
            else
                trd = squeeze(InputParam.Ts(:,1,:)); 
                cvd = squeeze(InputParam.Ts(:,2,:)); 
                tsd = squeeze(InputParam.Ts(:,3,:));
                if ~isempty(Yocv) 
                    if multoocv
                        ocv = [];
                        for n=1:numel(Yocv)
                            ocv = [ocv squeeze(InputParam.Ts(:,3+n,:))]; 
                        end
                    else
                        ocv = squeeze(InputParam.Ts(:,4,:));
                    end
                end
                if  any(any(cellfun(@isempty,trd))) ||  any(any(cellfun(@isempty,cvd))) ||  any(any(cellfun(@isempty,tsd)))
                    error('Empty training/test/validation data matrices returned by preprocessing pipeline. Check your settings')
                end
            end
            
            % Check whether we have imputed labels
            if isfield(SrcParam,'TrL_imputed')
                TrL = SrcParam.TrL_imputed; 
                [~,TrL] = nk_ManageNanCases(vTs{1}, TrL, SrcParam.iTr); 
                tTrL = labels(TrI,lb);
                TrL(~isnan(tTrL)) = tTrL(~isnan(tTrL));
            else
                % Overwrite labels adjusted to NaN cases
                TrL = labels(TrI(~SrcParam.iTr),lb); 
                if ~oocvonly
                    CVL = labels(CVI(~SrcParam.iCV),lb);
                end
            end
            
            clear InputParam
            switch BINMOD

                case 0 % Multi-group mode both in FBINMOD and PREPROC.BINMOD
                    if ukbin > 1
                        [tY.Tr{k,l}{u},TrL] = nk_ManageNanCases(trd, TrL, SrcParam.iTr);
                        if ~oocvonly
                            [tY.CV{k,l}{u},CVL] = nk_ManageNanCases(cvd, CVL, SrcParam.iCV);
                            tY.Ts{k,l}{u} = nk_ManageNanCases(tsd, [], SrcParam.iTs);
                        end
                        if ~isempty(Yocv) 
                            if multoocv
                                for n=1:numel(Yocv)
                                    tYocv.Ts{k,l}{u}{n} = nk_ManageNanCases(ocv{n},[],SrcParam.iOCV{n}); 
                                end
                            else
                                tYocv.Ts{k,l}{u} = nk_ManageNanCases(ocv,[],SrcParam.iOCV); 
                            end
                        end
                    else
                        tY.Tr{k,l} = nk_ManageNanCases(trd, TrL, SrcParam.iTr);
                        if ~oocvonly
                            [tY.CV{k,l},CVL] = nk_ManageNanCases(cvd, CVL, SrcParam.iCV);
                            tY.Ts{k,l} = nk_ManageNanCases(tsd, [], SrcParam.iTs);
                        end
                        if ~isempty(Yocv)
                            if multoocv
                                for n=1:numel(Yocv)
                                    tYocv.Ts{k,l}{n} = nk_ManageNanCases(ocv{n},[],SrcParam.iOCV{n}); 
                                end
                            else
                                tYocv.Ts{k,l} = nk_ManageNanCases(ocv,[],SrcParam.iOCV); 
                            end
                        end
                    end

                    oTrL = labels(TrInd(tCV.cvin{i,j}.TrainInd{k,l}),lb);
                    oCVL = labels(TrInd(tCV.cvin{i,j}.TestInd{k,l}),lb);

                    % Write dichotomization labels to CV1 partition
                    for zu=1:kbin % Binary loop depending on the # of binary comparisons
                        % Generate logical indices
                        if isfield(tCV,'class') && length(tCV.class{i,j}{zu}.groups) == 2
                            % One-vs-One
                            indtr = ( oTrL == tCV.class{i,j}{zu}.groups(1) | oTrL == tCV.class{i,j}{zu}.groups(2) ) ;
                            if ~oocvonly
                                indcv = ( oCVL == tCV.class{i,j}{zu}.groups(1) | oCVL == tCV.class{i,j}{zu}.groups(2) ) ;
                            end
                        else
                            % One-vs-All
                            indtr = oTrL~=0;
                            if ~oocvonly
                                indcv = oCVL~=0;
                            end
                        end
                        % Write indices
                        tY.TrInd{k,l}{zu} = indtr;
                        if ~oocvonly
                            tY.CVInd{k,l}{zu} = indcv;
                        end
                        % Write labels to CV1 partition
                        if isfield(tCV,'class')
                            tY.TrL{k,l}{zu} = tCV.class{i,j}{zu}.TrainLabel{k,l};
                            if ~oocvonly, tY.CVL{k,l}{zu} = tCV.class{i,j}{zu}.TestLabel{k,l}; end
                        else
                            tY.TrL{k,l}{zu} = labels(tCV.TrainInd{i,j}(tCV.cvin{i,j}.TrainInd{k,l}),lb);
                            if ~oocvonly, tY.CVL{k,l}{zu} = labels(tCV.TrainInd{i,j}(tCV.cvin{i,j}.TestInd{k,l}),lb); end
                        end
                    end

                case 1   

                    % Write data to CV1 partition
                    [tY.Tr{k,l}{u},TrL] = nk_ManageNanCases(trd, TrL, SrcParam.iTr); 
                    if ~oocvonly
                        [tY.CV{k,l}{u},CVL] = nk_ManageNanCases(cvd, CVL, SrcParam.iCV); 
                        tY.Ts{k,l}{u} = nk_ManageNanCases(tsd, [], SrcParam.iTs);
                    end
                    if ~isempty(Yocv)
                        if multoocv 
                            for n=1:numel(Yocv)
                                tYocv.Ts{k,l}{u}{n} = nk_ManageNanCases(ocv{n}, [], SrcParam.iOCV{n}); 
                            end
                        else
                            tYocv.Ts{k,l}{u} = nk_ManageNanCases(ocv, [], SrcParam.iOCV); 
                        end
                    end
                    
                    oTrL = labels(TrInd(tCV.cvin{i,j}.TrainInd{k,l}),lb);
                    oCVL = labels(TrInd(tCV.cvin{i,j}.TestInd{k,l}),lb);

                    if ~strcmp(MODEFL,'regression') && length(tCV.class{i,j}{u}.groups) == 2
                        indtr = ( oTrL == tCV.class{i,j}{u}.groups(1) | oTrL == tCV.class{i,j}{u}.groups(2) );
                        if ~oocvonly, indcv = ( oCVL == tCV.class{i,j}{u}.groups(1) | oCVL == tCV.class{i,j}{u}.groups(2) ) ; end
                    else
                        indtr = true(size(TrI,1),1);
                        if ~oocvonly, indcv = true(size(CVI,1),1); end
                    end
                    % Write indices
                    tY.TrInd{k,l}{u} = indtr;
                    if ~oocvonly, tY.CVInd{k,l}{u} = indcv; end

                    switch MODEFL
                        case 'regression' 
                            tY.TrL{k,l}{u} = oTrL;
                            if ~oocvonly, tY.CVL{k,l}{u} = oCVL; end
                        case 'classification'
                            if RAND.Decompose ~=9
                                % Write labels to CV1 partition
                                tY.TrL{k,l}{u} = tCV.class{i,j}{u}.TrainLabel{k,l}(:,lb);	
                                if ~oocvonly, tY.CVL{k,l}{u} = tCV.class{i,j}{u}.TestLabel{k,l}(:,lb); end
                            else
                                tY.TrL{k,l}{u} = TrL;
                                if ~oocvonly, tY.CVL{k,l}{u} = labels(tCV.TrainInd{i,j}(tCV.cvin{i,j}.TestInd{k,l}),lb); end
                            end
                    end

            end
            if ~isempty(MULTI) && MULTI.flag, tY.mTrL{k,l} = TrL; tY.mCVL{k,l} = CVL; end
            
            % save parameters
            if paramfl.write || cv2flag, Pnt(k,l,u).TrainedParam = oTrainedParam; end
            clear TrainedParam SrcParam
        end
        fprintf('\tCompleted in %1.2fs. ',toc(tElapsed)); 
    end
end

