function cv_ExportModels(inp)
% =========================================================================
% FORMAT Results = cv_ExportModels(inp)
% =========================================================================
% The export module saves optimally trained model and preprocessing structures 
% for each CV2 partition in a opt subdirectory of the given analysis
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris & Clara Vetter, last modified 08/2023
global SVM RFE MODEFL CV SCALE SAV CVPOS RAND

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
FullPartFlag    = RFE.ClassRetrain;

multiflag       = false; if strcmp(MODEFL,'classification'), multiflag = inp.multiflag; end
saveparam       = true; % inp.saveparam;
nclass          = inp.nclass;
analysis        = inp.analysis;
GridAct         = inp.GridAct;
batchflag       = inp.batchflag;
algostr         = GetMLType(SVM);


% Setup CV2 and CV1 counters and data containers:
[ix, jx]        = size(CV.TrainInd);
[iy, jy]        = size(CV.cvin{1,1}.TrainInd);
ll              = 1;
totLearn        = 0;

if ~exist('GridAct','var') || isempty(GridAct), GridAct = nk_CVGridSelector(ix,jx); end

% Train models with CV2 test data
if FullPartFlag
    % Use entire inner-cycle data to retrain models
    TrainWithCV2TsStr = 'CV1-Tr + CV1-Ts';
else
    % Use only inner-cycle training to retrain models
    TrainWithCV2TsStr = 'CV1-Tr';
end

% Check and transform labels if needed
inp = nk_ApplyLabelTransform( SCALE, MODEFL, inp );

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

% Parameter flag structure for preprocessing
paramfl = struct('use_exist',inp.loadparam, ...
    'found', false, ...
    'write', inp.saveparam, ...
    'writeCV1', inp.saveCV1, ...
    'multiflag', multiflag);

%Pre-smooth data, if needed, to save computational time
inp.ll=inp.GridAct';inp.ll=find(inp.ll(:));
inp = nk_PerfInitSpatial(analysis, inp, paramfl);

% =========================================================================
for f=1:ix % Loop through CV2 permutations
    for d=1:jx % Loop through CV2 folds
        fprintf('\n--------------------------------------------------------------------------')
        if ~GridAct(f,d)
            ll=ll+1;
            fprintf('\nSkipping CV2 [%g,%g] (user-defined).',f,d)
            continue
        end
        CVPOS.CV2p = f;
        CVPOS.CV2f = d;
        operm = f; ofold = d;
        inp.f = f; inp.d = d; inp.ll = ll;
        OptModelPath = nk_GenerateNMFilePath( inp.saveoptdir, SAV.matname, 'OptModel', [], inp.varstr, inp.id, operm, ofold);
        % Preprocess data at optimal hyperparameter combination(s) and save
        % trained hyperparameter structure
        [ inp, contfl, analysis, mapY, GD, MD, Param, paramfl, mapYocv ] = nk_ApplyTrainedPreproc(analysis, inp, paramfl);
        if contfl, continue; end
        fndMD = false;
        if ~fndMD, MD = cell(nclass,1); end

        % -----------------------------------------------------------------
        for h=1:nclass % Loop through binary comparisons

            if nclass > 1, fprintf('\n\n*** %s #%g ***',algostr, h); end

            %% Step 1: Get optimal model parameters
            % Retrieve optimal parameters from precomputed analysis structure
            % Differentiate according to binary or multi-group mode
            [Ps, Pspos, nP, Pdesc] = nk_GetModelParams2(analysis, multiflag, ll, h, inp.curlabel);

            if ~fndMD , MD{h} = cell(nP,1); end
            
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
                        [ TR , CV1, ~, ~ ] = nk_ReturnAtOptPos(mapY.Tr{k,l}{hix},  mapY.CV{k,l}{hix}, mapY.Ts{k,l}{hix}, mapYocv, Param{1}(k,l,hix), pnt);
                        if FullPartFlag, TR = [ TR; CV1]; end

                        % Get and build label info
                        modelTrL = mapY.TrL{k,l}{h};
                        if FullPartFlag, modelTrL = [modelTrL; mapY.CVL{k,l}{h}];    TrInd = [TrInd; CVInd]; end
                      
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

                        % Loop through feature subspaces
                        for u=1:ul
                            % Extract features according to mask
                            Ymodel = nk_ExtractFeatures(TR, F, [], u);
                            if ~fndMD
                                fprintf('Computing OptModel');
                                [~, MD{h}{m}{k,l}{u}] = nk_GetParam2(Ymodel, modelTrL, sPs, 1);
                            end
                        end
                    end
                end
            end
            % Now save trained model structure for current partition CV2 to disk
            if saveparam, fprintf('\nSaving %s', OptModelPath); save(OptModelPath, 'MD', 'ofold','operm'); end
        end
        ll=ll+1;
    end
end
