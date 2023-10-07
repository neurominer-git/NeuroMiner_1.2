function [InputParam, TrainedParam, SrcParam] = nk_GenPreprocSequence(InputParam, TemplParam, SrcParam, TrainedParam)
% =========================================================================
% FORMAT [InputParam, TrainedParam, SrcParam] = nk_GenPreprocSequence( ...
%                           InputParam, TemplParam, SrcParam, TrainedParam)
% =========================================================================
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 05/2022
global MODEFL NM VERBOSE

if isfield(TemplParam,'ACTPARAM')
    
    % Prepare preprocessing
    lact                    = numel(TemplParam.ACTPARAM);
    InputParam.P            = cell(lact,1);
    InputParam.CV1perm      = SrcParam.CV1perm;
    InputParam.CV1fold      = SrcParam.CV1fold;
    actionseq               = cell(1,lact);
    InputParam.curclass     = SrcParam.u;
    if ~isfield(SrcParam,'oocvonly'), SrcParam.oocvonly = 0; end
    if SrcParam.oocvonly, tscnt = 0; else, tscnt=3; end

    % Loop through ACTPARAM sequence
    for ac = 1:lact
        
        actionseq{ac} = TemplParam.ACTPARAM{ac}.cmd;
        InputParam.P{ac}.LabelInteraction = false;
        InputParam.P{ac}.cmd = actionseq{ac};
        
        switch actionseq{ac}
            
            case 'labelimpute'
                
                if VERBOSE, fprintf('\n* Impute NaNs in label vector.'); end
                InputParam.P{ac}.LABELIMPUTE = TemplParam.ACTPARAM{ac}.LABELIMPUTE;
                SrcParam.NaNflag = true;
                InputParam.P{ac}.BINMOD = TemplParam.BINMOD;

                    
            case 'impute'
                
                if VERBOSE, fprintf('\n* Impute NaNs in matrix block of %g columns.',sum(TemplParam.ACTPARAM{ac}.IMPUTE.blockind)); end
                InputParam.P{ac}.IMPUTE = TemplParam.ACTPARAM{ac}.IMPUTE;
        
            case 'copyfeat'

                if VERBOSE, fprintf('\n* Use all data.'); end
            
            case 'normalize'
                
                if ~isempty(SrcParam.covars)
                    InputParam.P{ac}.IND = TemplParam.ACTPARAM{ac}.IND;
                    InputParam.P{ac}.TsInd = [];
                    if VERBOSE
                        fprintf('\n* PER-GROUP NORMALIZATION AND ZERO-VARIANCE REMOVAL ACROSS GROUPS')
                        fprintf('\n\t- Normalize data to mean effects in group: %s', NM.covnames{TemplParam.ACTPARAM{ac}.IND});
                    end
                    InputParam.P{ac}.TrInd        = SrcParam.covars( SrcParam.TrX, TemplParam.ACTPARAM{ac}.IND );  
                    InputParam.P{ac}.TrInd(SrcParam.iTrX,:)=[]; 
                    if isfield(SrcParam,'TrI')         
                        InputParam.P{ac}.TsInd{end+1} = SrcParam.covars( SrcParam.TrI, TemplParam.ACTPARAM{ac}.IND );
                        InputParam.P{ac}.TsInd{1}(SrcParam.iTr,:)=[];
                    end
                    if isfield(SrcParam,'CVI')         
                        InputParam.P{ac}.TsInd{2} = SrcParam.covars( SrcParam.CVI, TemplParam.ACTPARAM{ac}.IND );
                        InputParam.P{ac}.TsInd{2}(SrcParam.iCV,:)=[];
                    end
                    if isfield(SrcParam,'TsI')         
                        InputParam.P{ac}.TsInd{3} = SrcParam.covars( SrcParam.TsI, TemplParam.ACTPARAM{ac}.IND );
                        InputParam.P{ac}.TsInd{3}(SrcParam.iTs,:)=[];
                    end
                    if ~isempty(SrcParam.covars_oocv) 
                        if iscell(SrcParam.covars_oocv)
                            for n=1:numel(SrcParam.covars_oocv)
                                InputParam.P{ac}.TsInd{tscnt+1} = SrcParam.covars_oocv{n}( :, TemplParam.ACTPARAM{ac}.IND );   
                                InputParam.P{ac}.TsInd{tscnt+n}(SrcParam.iOCV{n},:)=[]; 
                            end
                        else
                            InputParam.P{ac}.TsInd{tscnt+1} = SrcParam.covars_oocv( :, TemplParam.ACTPARAM{ac}.IND ); 
                            InputParam.P{ac}.TsInd{tscnt+1}(SrcParam.iOCV,:)=[]; 
                        end
                    end
                end
                
            case 'correctnuis'
                
                if ~isempty(SrcParam.covars)
                    InputParam.P{ac}.COVAR = TemplParam.ACTPARAM{ac}.COVAR;
                    if VERBOSE 
                        fprintf('\n* ADJUSTING DATA FOR COVARIATE EFFECTS')
                        fprintf('\n\t- Nuisance covariate: %s', NM.covnames{TemplParam.ACTPARAM{ac}.COVAR})
                    end
                    InputParam.P{ac}.TsCovars = [];
                    if isfield(SrcParam,'TrX')         
                        InputParam.P{ac}.TrCovars        = SrcParam.covars( SrcParam.TrX, TemplParam.ACTPARAM{ac}.COVAR );
                        InputParam.P{ac}.TrCovars(SrcParam.iTrX,:)=[];
                    end
                    if isfield(SrcParam,'TrI')         
                        InputParam.P{ac}.TsCovars{end+1} = SrcParam.covars( SrcParam.TrI, TemplParam.ACTPARAM{ac}.COVAR );
                        InputParam.P{ac}.TsCovars{1}(SrcParam.iTr,:) = [];
                    end
                    if isfield(SrcParam,'CVI')         
                        InputParam.P{ac}.TsCovars{end+1} = SrcParam.covars( SrcParam.CVI, TemplParam.ACTPARAM{ac}.COVAR );   
                        InputParam.P{ac}.TsCovars{2}(SrcParam.iCV,:)=[];
                    end
                    if isfield(SrcParam,'TsI')         
                        InputParam.P{ac}.TsCovars{end+1} = SrcParam.covars( SrcParam.TsI, TemplParam.ACTPARAM{ac}.COVAR );
                        InputParam.P{ac}.TsCovars{3}(SrcParam.iTs,:)=[]; 
                    end
                    if ~isempty(SrcParam.covars_oocv)  
                        if iscell(SrcParam.covars_oocv)
                            for n=1:numel(SrcParam.covars_oocv)
                                InputParam.P{ac}.TsCovars{tscnt+n} = SrcParam.covars_oocv{n}( :, TemplParam.ACTPARAM{ac}.COVAR );
                                if iscell(SrcParam.iOCV)
                                    InputParam.P{ac}.TsCovars{tscnt+n}(SrcParam.iOCV{n},:)=[]; 
                                else
                                    InputParam.P{ac}.TsCovars{tscnt+n}(SrcParam.iOCV,:)=[]; 
                                end
                            end
                        else
                            InputParam.P{ac}.TsCovars{tscnt+1} = SrcParam.covars_oocv( :, TemplParam.ACTPARAM{ac}.COVAR );
                            InputParam.P{ac}.TsCovars{tscnt+1}(SrcParam.iOCV,:)=[]; 
                        end
                    end
                    if isfield(TemplParam.ACTPARAM{ac},'METHOD') && TemplParam.ACTPARAM{ac}.METHOD==2
                            if VERBOSE, fprintf('\n\t- Method: Combat'); end
                            InputParam.P{ac}.METHOD = 2; InputParam.P{ac}.TsMod = [];
                            if TemplParam.ACTPARAM{ac}.MCOVARLABEL
                                covars = [ SrcParam.covars(:, TemplParam.ACTPARAM{ac}.MCOVAR) NM.label ];
                                if ~isempty(SrcParam.covars_oocv) 
                                    if iscell(SrcParam.covars_oocv)
                                        for n=1:numel(SrcParam.covars_oocv)
                                            covars_oocv{n} = [ SrcParam.covars_oocv{n}(:,TemplParam.ACTPARAM{ac}.MCOVAR) ...
                                                zeros(size(SrcParam.covars_oocv{n},1), size(dummy_labels,2)) ];
                                        end
                                    else
                                        covars_oocv = [ SrcParam.covars_oocv(:,TemplParam.ACTPARAM{ac}.MCOVAR) ...
                                            zeros(size(SrcParam.covars_oocv,1), size(dummy_labels,2)) ];
                                    end
                                end
                            else
                                covars = SrcParam.covars( : , TemplParam.ACTPARAM{ac}.MCOVAR ); 
                                if ~isempty(SrcParam.covars_oocv)
                                    if iscell(SrcParam.covars_oocv)
                                        for n=1:numel(SrcParam.covars_oocv)
                                            covars_oocv{n} = SrcParam.covars_oocv{n}( : , TemplParam.ACTPARAM{ac}.MCOVAR );
                                        end
                                    else
                                        covars_oocv = SrcParam.covars_oocv ( : , TemplParam.ACTPARAM{ac}.MCOVAR ); 
                                    end
                                end
                               
                            end
                            InputParam.P{ac}.TrMod        = covars( SrcParam.TrX, :);
                            InputParam.P{ac}.TrMod(SrcParam.iTrX,:)   = []; 
                            
                            if isfield(SrcParam,'TrI')         
                                InputParam.P{ac}.TsMod{1} = covars( SrcParam.TrI, :);
                                InputParam.P{ac}.TsMod{1}(SrcParam.iTr,:) = [];
                            end
                            if isfield(SrcParam,'CVI')         
                                InputParam.P{ac}.TsMod{2} = covars( SrcParam.CVI, :);
                                InputParam.P{ac}.TsMod{2}(SrcParam.iCV,:)=[]; 
                            end
                            if isfield(SrcParam,'TsI')         
                                InputParam.P{ac}.TsMod{3} = covars( SrcParam.TsI, :); 
                                InputParam.P{ac}.TsMod{3}(SrcParam.iTs,:)=[];   
                            end
                            if ~isempty(SrcParam.covars_oocv) 
                                if iscell(covars_oocv)
                                    for n=1:numel(SrcParam.covars_oocv)
                                        InputParam.P{ac}.TsMod{tscnt+n} = covars_oocv{n};
                                        InputParam.P{ac}.TsMod{tscnt+n}(SrcParam.iOCV{n},:)=[]; 
                                    end
                                else
                                    InputParam.P{ac}.TsMod{tscnt+1} = covars_oocv;
                                    InputParam.P{ac}.TsMod{tscnt+1}(SrcParam.iOCV,:)=[]; 
                                end
                            end
                    elseif isfield(TemplParam.ACTPARAM{ac},'METHOD') && TemplParam.ACTPARAM{ac}.METHOD==3
                        if VERBOSE, fprintf('\n\t- Method: fastICA analysis'); end
                        InputParam.P{ac}.METHOD = 3;
                    else
                        if VERBOSE, fprintf('\n\t- Method: Partial correlations analysis'); end
                        InputParam.P{ac}.METHOD = 1;
                    end
                    
                end
                if InputParam.P{ac}.METHOD == 1
                    if isfield(TemplParam.ACTPARAM{ac},'INTERCEPT')
                        InputParam.P{ac}.INTERCEPT = TemplParam.ACTPARAM{ac}.INTERCEPT-1;
                        if VERBOSE
                            switch InputParam.P{ac}.INTERCEPT
                                case 0
                                    fprintf('\n\t-> Not including intercept.'); 
                                case 1
                                    fprintf('\n\t-> Including intercept.'); 
                            end  
                        end
                    end
                    if isfield(TemplParam.ACTPARAM{ac},'COVDIR')
                        InputParam.P{ac}.COVDIR = TemplParam.ACTPARAM{ac}.COVDIR-1;
                        if VERBOSE 
                            switch InputParam.P{ac}.COVDIR
                                case 0
                                    fprintf('\n\t-> Covariate effects will be removed from data.'); 
                                case 1
                                    fprintf('\n\t-> Covariate effects will be increased in data.')
                            end  
                        end
                    end
                    if isfield(TemplParam.ACTPARAM{ac},'BETAEXT') && ~isempty(TemplParam.ACTPARAM{ac}.BETAEXT)
                        InputParam.P{ac}.BETAEXT = TemplParam.ACTPARAM{ac}.BETAEXT;
                        if VERBOSE,fprintf('\n\t Beta parameter(s) computed in an OOT-sample will be used.'); end
                    else
                        if isfield(TemplParam.ACTPARAM{ac},'SUBGROUP') && ~isempty(TemplParam.ACTPARAM{ac}.SUBGROUP) && ~isnan(any(TemplParam.ACTPARAM{ac}.SUBGROUP))
                            InputParam.P{ac}.SUBGROUP = TemplParam.ACTPARAM{ac}.SUBGROUP(SrcParam.TrX,:);
                            InputParam.P{ac}.SUBGROUP(SrcParam.iTrX) = []; 
                            if VERBOSE,fprintf('\n\t-> Beta parameter(s) will be computed from a specific subgroup.'); end
                        end
                    end
                elseif InputParam.P{ac}.METHOD == 2
                    if isfield(TemplParam.ACTPARAM{ac},'SUBGROUP') && ~isempty(TemplParam.ACTPARAM{ac}.SUBGROUP)
                        InputParam.P{ac}.SUBGROUP = TemplParam.ACTPARAM{ac}.SUBGROUP(SrcParam.TrX,:);
                        InputParam.P{ac}.SUBGROUP(SrcParam.iTrX) = [];
                        if VERBOSE,fprintf('\n\t-> Combat parameter(s) will be computed from a specific subgroup.'); end
                    end
                     InputParam.P{ac}.COVDIR=0;
                     InputParam.P{ac}.INTERCEPT=0;
                elseif InputParam.P{ac}.METHOD == 3
                    if isfield(TemplParam.ACTPARAM{ac},'SUBGROUP') && ~isempty(TemplParam.ACTPARAM{ac}.SUBGROUP)
                        InputParam.P{ac}.SUBGROUP = TemplParam.ACTPARAM{ac}.SUBGROUP(SrcParam.TrX,:);
                        InputParam.P{ac}.SUBGROUP(SrcParam.iTrX) = [];
                        if VERBOSE,fprintf('\n\t-> FastICA parameter(s) will be computed from a specific subgroup.'); end
                    end
                     InputParam.P{ac}.COVDIR=0;
                     InputParam.P{ac}.INTERCEPT=0;
                end
                if isfield(TemplParam.ACTPARAM{ac},'featind')
                     if VERBOSE, fprintf('\n\t- Feature subspace for covariate correction identified.'); end
                     InputParam.P{ac}.featind = TemplParam.ACTPARAM{ac}.featind;
                end
                
            case 'remmeandiff'
                
                if isfield(NM,'covars') && ~isempty(NM.covars)
                    InputParam.P{ac}.sIND = TemplParam.ACTPARAM{ac}.sIND; InputParam.P{ac}.dIND = TemplParam.ACTPARAM{ac}.dIND;
                    if VERBOSE
                        fprintf('\n* OFFSET CORRECTION')
                        fprintf('\n\t - Compute group mean offset of : %s', NM.covnames{TemplParam.ACTPARAM{ac}.sIND})
                        fprintf('\n\t - Remove global mean and offsets from : %s', NM.covnames{TemplParam.ACTPARAM{ac}.dIND})
                    end
                    if isfield(SrcParam,'TrX')         
                        InputParam.P{ac}.sTrInd        = SrcParam.covars( SrcParam.TrX, TemplParam.ACTPARAM{ac}.sIND );
                        InputParam.P{ac}.dTrInd        = SrcParam.covars( SrcParam.TrX, TemplParam.ACTPARAM{ac}.dIND );
                        InputParam.P{ac}.sTrInd( SrcParam.iTrX,:) = []; 
                        InputParam.P{ac}.dTrInd( SrcParam.iTrX,:) = []; 
                    end
                    if isfield(SrcParam,'TrI')         
                        InputParam.P{ac}.sTsInd{1}     = SrcParam.covars( SrcParam.TrI, TemplParam.ACTPARAM{ac}.sIND );   
                        InputParam.P{ac}.dTsInd{1}     = SrcParam.covars( SrcParam.TrI, TemplParam.ACTPARAM{ac}.dIND );
                        InputParam.P{ac}.sTsInd{1}( SrcParam.iTr,:) = []; 
                        InputParam.P{ac}.dTsInd{1}( SrcParam.iTr,:) = []; 
                    end
                    if isfield(SrcParam,'CVI')         
                        InputParam.P{ac}.sTsInd{end+1} = SrcParam.covars( SrcParam.CVI, TemplParam.ACTPARAM{ac}.sIND );
                        InputParam.P{ac}.dTsInd{end+1} = SrcParam.covars( SrcParam.CVI, TemplParam.ACTPARAM{ac}.dIND );
                        InputParam.P{ac}.sTsInd{2}( SrcParam.iCV,:) = []; 
                        InputParam.P{ac}.dTsInd{2}( SrcParam.iCV,:) = []; 
                    end
                    if isfield(SrcParam,'TsI')         
                        InputParam.P{ac}.sTsInd{end+1} = SrcParam.covars( SrcParam.TsI, TemplParam.ACTPARAM{ac}.sIND );
                        InputParam.P{ac}.dTsInd{end+1} = SrcParam.covars( SrcParam.TsI, TemplParam.ACTPARAM{ac}.dIND );
                        InputParam.P{ac}.sTsInd{3}( SrcParam.iTs,:) = []; 
                        InputParam.P{ac}.dTsInd{3}( SrcParam.iTs,:) = []; 
                        
                    end
                    if ~isempty(SrcParam.covars_oocv)
                        if iscell(SrcParam.covars)
                            for n=1:numel(SrcParam.covars_oocv)
                                InputParam.P{ac}.sTsInd{tscnt+n} = SrcParam.covars_oocv{n}( : , TemplParam.ACTPARAM{ac}.sIND );   
                                InputParam.P{ac}.dTsInd{tscnt+n} = SrcParam.covars_oocv{n}( : , TemplParam.ACTPARAM{ac}.dIND );
                                InputParam.P{ac}.sTsInd{tscnt+n}(SrcParam.iOCV{n},:) = [];  
                                InputParam.P{ac}.dTsInd{tscnt+n}(SrcParam.iOCV{n},:) = [];
                            end
                        else
                            InputParam.P{ac}.sTsInd{tscnt+1} = SrcParam.covars_oocv( : , TemplParam.ACTPARAM{ac}.sIND );   
                            InputParam.P{ac}.dTsInd{tscnt+1} = SrcParam.covars_oocv( : , TemplParam.ACTPARAM{ac}.dIND );
                            InputParam.P{ac}.sTsInd{tscnt+1}(SrcParam.iOCV,:) = [];  
                            InputParam.P{ac}.dTsInd{tscnt+1}(SrcParam.iOCV,:) = [];
                        end
                    end
                end
                
            case 'standardize' 
                
                % **************** STANDARDIZATION ***************
               
                if VERBOSE, fprintf('\n* Z-NORMALIZATION'); end
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX.opt)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                end
                if isfield(TemplParam.ACTPARAM{ac},'METHOD') && ~isempty(TemplParam.ACTPARAM{ac}.METHOD)
                    InputParam.P{ac}.method = TemplParam.ACTPARAM{ac}.METHOD;
                end
                if isfield(TemplParam.ACTPARAM{ac},'sIND') && ~isempty(TemplParam.ACTPARAM{ac}.sIND)
                    InputParam.P{ac}.sIND = TemplParam.ACTPARAM{ac}.sIND ;
                    if isfield(TemplParam.ACTPARAM{ac},'CALIBUSE') && TemplParam.ACTPARAM{ac}.CALIBUSE
                        InputParam.P{ac}.CALIBUSE      = true;
                    else
                        InputParam.P{ac}.CALIBUSE      = false;
                    end
                    if isfield(SrcParam,'TrX')  
                        InputParam.P{ac}.sTrInd        = SrcParam.covars( SrcParam.TrX, TemplParam.ACTPARAM{ac}.sIND );
                        InputParam.P{ac}.sTrInd(SrcParam.iTrX,:)=[];
                        denom = numel(InputParam.P{ac}.sTrInd);
                    elseif isfield(SrcParam,'TrI')
                        InputParam.P{ac}.sTsInd{1}     = SrcParam.covars( SrcParam.TrI, TemplParam.ACTPARAM{ac}.sIND );
                        InputParam.P{ac}.sTsInd{1}(SrcParam.iTr,:)=[]; 
                        denom = numel(InputParam.P{ac}.sTsInd{1});
                    end
                     if VERBOSE,fprintf('\n\tMean & SD computation restricted to N=%g (%1.0f%% of training sample).', ...
                        sum(InputParam.P{ac}.sTrInd), sum(InputParam.P{ac}.sTrInd) * 100 / denom); end
                else
                    InputParam.P{ac}.sTrInd = [];
                end
                
                %InputParam.P{ac}.dTsInd = cell(4,1);
                if isfield(TemplParam.ACTPARAM{ac},'dIND') && ~isempty(TemplParam.ACTPARAM{ac}.dIND)  
                    InputParam.P{ac}.dIND = TemplParam.ACTPARAM{ac}.dIND ;
                    if isfield(SrcParam,'TrX')         
                        InputParam.P{ac}.dTrInd    = SrcParam.covars( SrcParam.TrX, TemplParam.ACTPARAM{ac}.dIND );
                        InputParam.P{ac}.dTrInd(SrcParam.iTrX,:)=[];
                    end
                    if isfield(SrcParam,'TrI')         
                        InputParam.P{ac}.dTsInd{1} = SrcParam.covars( SrcParam.TrI, TemplParam.ACTPARAM{ac}.dIND );
                        InputParam.P{ac}.dTsInd{1}(SrcParam.iTr,:)=[];
                    end
                    if isfield(SrcParam,'CVI')         
                        InputParam.P{ac}.dTsInd{2} = SrcParam.covars( SrcParam.CVI, TemplParam.ACTPARAM{ac}.dIND );
                        InputParam.P{ac}.dTsInd{2}(SrcParam.iCV,:)=[];
                    end
                    if isfield(SrcParam,'TsI')         
                        InputParam.P{ac}.dTsInd{3} = SrcParam.covars( SrcParam.TsI, TemplParam.ACTPARAM{ac}.dIND );
                        InputParam.P{ac}.dTsInd{3}(SrcParam.iTs,:)=[]; 
                    end
                    if ~isempty(SrcParam.covars_oocv) 
                        if iscell(SrcParam.covars_oocv)
                            for n=1:numel(SrcParam.covars_oocv)
                                InputParam.P{ac}.dTsInd{tscnt+n} = SrcParam.covars_oocv{n}(:,TemplParam.ACTPARAM{ac}.dIND);
                                InputParam.P{ac}.dTsInd{tscnt+n}(SrcParam.iOCV{n},:)=[]; 
                            end
                        else
                            InputParam.P{ac}.dTsInd{tscnt+1} = SrcParam.covars_oocv(:, TemplParam.ACTPARAM{ac}.dIND);
                            InputParam.P{ac}.dTsInd{tscnt+1}(SrcParam.iOCV,:)=[]; 
                        end
                    end
                end
                
                InputParam.P{ac}.WINSOPT = TemplParam.ACTPARAM{ac}.WINSOPT;
                if isfield(TemplParam.ACTPARAM{ac},'IQRFUN')
                    InputParam.P{ac}.IQRFUN = TemplParam.ACTPARAM{ac}.IQRFUN;
                end
                
            case 'spatialfilter'
                
                if VERBOSE, fprintf('\n* SPATIAL FILTERING'); end
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX.opt)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                end
                InputParam.P{ac}.SPATIAL = TemplParam.ACTPARAM{ac}.SPATIAL;
                
            case 'unitnormalize'
                % **************** UNIT VECTOR NORMALIZATION ***************
                if VERBOSE, fprintf('\n* UNIT VECTOR NORMALIZATION'); end
                InputParam.P{ac}.METHOD = TemplParam.ACTPARAM{ac}.METHOD;
                
            case 'discretize'
                
                % **************** DISCRETIZATION *****************
                if VERBOSE,fprintf('\n* DISCRETIZATION'); end                   
                InputParam.P{ac}.DISCRET = TemplParam.ACTPARAM{ac}.DISCRET;
                InputParam.P{ac}.method = str2func('PerfDiscretizeObj');
                
            case 'symbolize'
                
                % **************** SYMBOLIZATION *****************
                if VERBOSE,fprintf('\n* SYMBOLIZATION'); end               
                InputParam.P{ac}.SYMBOL = TemplParam.ACTPARAM{ac}.SYMBOL;
                InputParam.P{ac}.method = str2func('PerfSymbolizeObj');
                                
            case 'reducedim'
                
                % **************** DIM. REDUCTION *****************
                if ~strcmp(TemplParam.ACTPARAM{ac}.DR.RedMode,'none')
                    if VERBOSE,fprintf('\n* DIMENSIONALITY REDUCTION'); end
                    InputParam.P{ac}.DR = TemplParam.ACTPARAM{ac}.DR;
                    InputParam.P{ac}.DR.Modus = 'JDQR';
                    switch MODEFL
                        case 'classification' 
                            if TemplParam.BINMOD
                                InputParam.P{ac}.DR.labels = SrcParam.BinaryTrainLabel;
                            else
                                InputParam.P{ac}.DR.labels = SrcParam.MultiTrainLabel;
                            end
                        case 'regression'
                            InputParam.P{ac}.DR.labels = SrcParam.TrainLabel;
                    end
                    if strcmp(TemplParam.ACTPARAM{ac}.DR.RedMode,'PLS')
                        InputParam.P{ac}.DR.PLS.V = TemplParam.ACTPARAM{ac}.DR.PLS.V(SrcParam.TrX,:);
                        %InputParam.P{ac}.DR.PLS.VT{1} = TemplParam.ACTPARAM{ac}.DR.PLS.V(SrcParam.TrI,:);
                        InputParam.P{ac}.DR.PLS.VT{1} = TemplParam.ACTPARAM{ac}.DR.PLS.V(SrcParam.CVI,:);
                        InputParam.P{ac}.DR.PLS.VT{2} = TemplParam.ACTPARAM{ac}.DR.PLS.V(SrcParam.TsI,:);
                        if sum(SrcParam.iTrX), InputParam.P{ac}.DR.PLS.V(SrcParam.iTrX,:)=[]; end
                        %if sum(SrcParam.iTr), InputParam.P{ac}.DR.PLS.VT{1}(SrcParam.iTr,:)=[]; end
                        if sum(SrcParam.iCV), InputParam.P{ac}.DR.PLS.VT{1}(SrcParam.iCV,:)=[]; end
                        if sum(SrcParam.iTs), InputParam.P{ac}.DR.PLS.VT{2}(SrcParam.iTs,:)=[]; end
                    end
                    if strcmp(TemplParam.ACTPARAM{ac}.DR.RedMode,{'PLS','LDA','KLDA', 'KFDA', 'KernelLDA', 'KernelFDA', 'GDA', 'NCA', 'LMNN'})
                        InputParam.P{ac}.LabelInteraction = true;
                    end
                else
                    InputParam.P{ac}.DR.RedMode = 'none';
                    InputParam.P{ac}.DR.dims = size(InputParam.Tr,2);
                end
                
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                    PX = nk_ReturnParamChain(TemplParam.ACTPARAM{ac});
                    InputParam.P{ac}.DR.Params = PX.Params;
                    InputParam.P{ac}.DR.Params_desc = PX.Params_desc;
                end
                
            case 'scale'
                
                % ******************* SCALING *********************
                if isfield(TemplParam.ACTPARAM{ac}.SCALE,'ZeroOne')
                    if VERBOSE
                        switch TemplParam.ACTPARAM{ac}.SCALE.ZeroOne
                            case 1
                                fprintf('\n* SCALING [0 1]')
                            case 2
                                fprintf('\n* SCALING [-1 1]')
                        end
                    end
                end
                InputParam.P{ac}.SCALE = TemplParam.ACTPARAM{ac}.SCALE;
           
            case 'rankfeat'
                
                % **************** FEATURE RANKING ***************
                if VERBOSE, fprintf('\n* FEATURE WEIGHTING'); end
                InputParam.P{ac}.RANK = TemplParam.ACTPARAM{ac}.RANK;
                InputParam.P{ac}.RANK.curlabel = TemplParam.ACTPARAM{ac}.RANK.label(SrcParam.TrX,:);
                InputParam.P{ac}.RANK.curlabel(SrcParam.iTrX,:)=[]; 
                if isfield( TemplParam.ACTPARAM{ac}.RANK,'glabel' )
                    % glabel is a logical vector
                    InputParam.P{ac}.RANK.curglabel = TemplParam.ACTPARAM{ac}.RANK.glabel(SrcParam.TrX, :);
                    InputParam.P{ac}.RANK.curglabel(SrcParam.iTrX,:)=[]; 
                end
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                    PX = nk_ReturnParamChain(TemplParam.ACTPARAM{ac});
                    InputParam.P{ac}.RANK.Params = PX.Params;
                    InputParam.P{ac}.RANK.Params_desc = PX.Params_desc;
                end
                if strcmp(TemplParam.ACTPARAM{ac}.RANK.labeldesc,'NM target label') && ...
                        ~strcmp(TemplParam.ACTPARAM{ac}.RANK.algostr,'extern') && ...
                        ~strcmp(TemplParam.ACTPARAM{ac}.RANK.algostr,'extern2') 
                    InputParam.P{ac}.LabelInteraction = true;
                end
                
            case 'extfeat'
                % **************** FEATURE EXTRACTION ***************
                if VERBOSE, fprintf('\n* WEIGHTING-BASED FEATURE GENERATION'); end
                InputParam.P{ac}.W_ACT = TemplParam.ACTPARAM{ac}.W_ACT;
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                    PX = nk_ReturnParamChain(TemplParam.ACTPARAM{ac});
                    InputParam.P{ac}.W_ACT.Params = PX.Params;
                    InputParam.P{ac}.W_ACT.Params_desc = PX.Params_desc;
                end
                
            case 'extdim'
                 % ******************* COMPONENT EXTRACTION (HELPER OF DIM. REDUCTION) *********************
                if VERBOSE, fprintf('\n* COMPONENT EXTRACTION FOLLOWING DIMENSIONALITY REDUCTION'); end
                InputParam.P{ac}.EXTDIM = TemplParam.ACTPARAM{ac}.EXTDIM;
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                    PX = nk_ReturnParamChain(TemplParam.ACTPARAM{ac});
                    InputParam.P{ac}.EXTDIM.Params = PX.Params;
                    InputParam.P{ac}.EXTDIM.Params_desc = PX.Params_desc;
                end
                
            case 'elimzero'
                 % ******************* ELIMINATE ZERO VAR, NAN/INF FEATURES *********************
                if VERBOSE, fprintf('\n* PRUNE ATTRIBUTES'); end
                InputParam.P{ac}.PRUNE =  TemplParam.ACTPARAM{ac}.PRUNE;
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX.opt)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                end
                
            case 'remvarcomp'
                 % ******************* REMOVE VARIANCE COMPONENTS  *********************
                if VERBOSE, fprintf('\n* EXTRACT VARIANCE COMPONENTS'); end
                InputParam.P{ac}.REMVARCOMP = TemplParam.ACTPARAM{ac}.REMVARCOMP;
                InputParam.P{ac}.REMVARCOMP.G = TemplParam.ACTPARAM{ac}.REMVARCOMP.G(SrcParam.TrX,:);
                if ~isempty(SrcParam.iTrX) 
                    InputParam.P{ac}.REMVARCOMP.G(SrcParam.iTrX,:)=[]; 
                end
                if isfield(TemplParam.ACTPARAM{ac}.REMVARCOMP,'SUBGROUP') 
                    switch TemplParam.ACTPARAM{ac}.REMVARCOMP.SUBGROUP.flag
                        case 1
                            InputParam.P{ac}.REMVARCOMP.indX = true(numel(find(SrcParam.TrX)),1);
                        case 2
                            InputParam.P{ac}.REMVARCOMP.indX = TemplParam.ACTPARAM{ac}.REMVARCOMP.SUBGROUP.ind(SrcParam.TrX);
                        case 3
                            n = round(numel(find(SrcParam.TrX))/100*TemplParam.ACTPARAM{ac}.REMVARCOMP.SUBGROUP.indperc);
                            ind = true(numel(find(SrcParam.TrX)),1); randind = randperm(numel(ind),n);
                            InputParam.P{ac}.REMVARCOMP.indX = ind(randind);
                        case 4
                            InputParam.P{ac}.REMVARCOMP.indX = true(numel(InputParam.C,1));
                        case 5
                            InputParam.P{ac}.REMVARCOMP.indX = TemplParam.ACTPARAM{ac}.REMVARCOMP.SUBGROUP.ind;
                        case 6
                            n = size(InputParam.C,1)/100*TemplParam.ACTPARAM{ac}.REMVARCOMP.SUBGROUP.indperc;
                            ind = true(size(InputParam.C,1),1); randind = randperm(numel(ind),n);
                            InputParam.P{ac}.REMVARCOMP.indX = ind(randind);
                    end
                    InputParam.P{ac}.REMVARCOMP.indX(SrcParam.iTrX)=[]; 
                end
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                    PX = nk_ReturnParamChain(TemplParam.ACTPARAM{ac});
                    InputParam.P{ac}.REMVARCOMP.Params = PX.Params;
                    InputParam.P{ac}.REMVARCOMP.Params_desc = PX.Params_desc;
                end
                
            case 'devmap'
                 % ******************* SENSITIZE FEATURES TO DEVIATION FROM NORMATIVE MODEL *********************
                InputParam.P{ac}.DEVMAP = TemplParam.ACTPARAM{ac}.DEVMAP;
              
                if isfield(SrcParam,'TrX')         
                    InputParam.P{ac}.TrInd    = SrcParam.TrX;
                    InputParam.P{ac}.TrInd(SrcParam.iTrX,:)=[]; 
                end
                if isfield(SrcParam,'TrI')     
                    InputParam.P{ac}.TsInd{1} = SrcParam.TrI ;
                    InputParam.P{ac}.TsInd{1}(SrcParam.iTr,:)=[]; 
                end
                if isfield(SrcParam,'CVI')         
                    InputParam.P{ac}.TsInd{2} = SrcParam.CVI;
                    InputParam.P{ac}.TsInd{2}(SrcParam.iCV,:)=[]; 
                end
                if isfield(SrcParam,'TsI')         
                    InputParam.P{ac}.TsInd{3} = SrcParam.TsI;
                    InputParam.P{ac}.TsInd{3}(SrcParam.iTs,:)=[]; 
                end
                if numel(InputParam.Ts)>3
                    for n=1:numel(InputParam.Ts)-tscnt
                        InputParam.P{ac}.TsInd{tscnt+n} = true(size(InputParam.Ts{tscnt+n},1),1);
                    end
                end
                
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                    PX = nk_ReturnParamChain(TemplParam.ACTPARAM{ac});
                    InputParam.P{ac}.DEVMAP.Params = PX.Params;
                    InputParam.P{ac}.DEVMAP.Params_desc = PX.Params_desc;
                end
                
            case 'graphSparsity'
                if VERBOSE, fprintf('\n* APPLY SPARSITY THRESHOLD TO CONNECTIVITY MATRICES'); end
                InputParam.P{ac} =  TemplParam.ACTPARAM{ac};
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX.opt)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                end
             case 'graphMetrics'
                if VERBOSE, fprintf('\n* COMPUTE GRAPH METRICS FROM CONNECTIVITY MATRICES'); end
                InputParam.P{ac} =  TemplParam.ACTPARAM{ac};
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX.opt)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                end
            case 'graphComputation'
                if VERBOSE, fprintf('\n* COMPUTE CONNECTIVITY MATRICES OF INDIVIDUAL NETWORKS'); end
                InputParam.P{ac} =  TemplParam.ACTPARAM{ac};
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX.opt)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                end
            case 'customPreproc'
                if VERBOSE, fprintf('\n* APPLY CUSTOM PREPROCESSING STEP'); end
                InputParam.P{ac} =  TemplParam.ACTPARAM{ac};
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX.opt)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                end
            case 'JuSpace'
                if VERBOSE, fprintf('\n* JUSPACE COMPUTATION'); end
                InputParam.P{ac} =  TemplParam.ACTPARAM{ac};
                InputParam.P{ac}.JUSPACE.analdir = NM.analysis{1,NM.runtime.curanal}.parentdir;
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX.opt)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                end
            case 'ROImeans'
                if VERBOSE, fprintf('\n* ROI MEANS COMPUTATION'); end
                InputParam.P{ac} =  TemplParam.ACTPARAM{ac};
                if isfield(TemplParam.ACTPARAM{ac},'PX') && ~isempty(TemplParam.ACTPARAM{ac}.PX.opt)
                    InputParam.P{ac}.opt = TemplParam.ACTPARAM{ac}.PX.opt;
                end
       end
    end
    if VERBOSE, fprintf('\nPreprocessing sequence setup completed. Executing ...'); end
    if exist('TrainedParam','var') && ~isempty(TrainedParam)
        [SrcParam, InputParam, TrainedParam] = nk_PerfPreprocessObj_core(SrcParam, InputParam, TrainedParam, actionseq);
    else
        [SrcParam, InputParam, TrainedParam] = nk_PerfPreprocessObj_core(SrcParam, InputParam, [], actionseq);
    end
    clear actionseq
else
    if VERBOSE,fprintf('\nNo preprocessing sequence detected. Aborting ...'); end
    TrainedParam = [];
end

end
