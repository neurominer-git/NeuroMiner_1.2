function [act, varind] = nk_TrainClass_config(act, varind, parentstr)
% =========================================================================
% FORMAT [act, varind] = nk_TrainClass_config(act, varind, parentstr)
% =========================================================================
% Interface for defining the parameters of the training and prediction
% process. The function also creates default parameters if a new modality
% has been added to the NM workspace
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 08/2023

global NM CALIBAVAIL

if ~exist('act','var'), act = []; end
menustr = []; menuact = [];

%% Select variate to work on
if ~exist('varind','var') || isempty(varind)
    if isfield(NM,'TrainParam')
        varind = NM.TrainParam.FUSION.M(1);
    else
        varind = 1;
    end
end

if isfield(NM, 'TrainParam') && isfield(NM.TrainParam, 'LABEL') && NM.TrainParam.LABEL.flag
    modeflag = NM.TrainParam.LABEL.newmode;
else
    modeflag = NM.modeflag;
end

%% Check whether TrainParams already exist
if ~isfield(NM,'TrainParam')
    % Create default NM parameters space
    nk_CVpartition_config(true);
    NM.TrainParam.STACKING.flag = 2;
    NM.TrainParam.FUSION.flag   = 0;
    NM.TrainParam.FUSION.M      = 1;
    [~,NM.TrainParam.SVM]       = nk_LIBSVM_config(NM,[],1);
    NM.TrainParam.SVM.prog      = 'LIBSVM';
    NM.TrainParam.SVM           = nk_Kernel_config(NM.TrainParam.SVM,1);
    NM.TrainParam.SVM.GridParam = 1;
    if strcmp(modeflag, 'regression'), NM.TrainParam.SVM.GridParam = 18; end
    NM.TrainParam.MULTI.flag    = 0;
    NM.TrainParam               = nk_Grid_config(NM.TrainParam, NM.TrainParam.SVM, varind, true);
    [~,NM.TrainParam.RFE]       = nk_RFE_config([], NM.TrainParam, NM.TrainParam.SVM, modeflag, NM.TrainParam.MULTI, NM.TrainParam.GRD, 1);
    NM.TrainParam.verbosity     = 1;
    NM.TrainParam.LABEL.flag    = 0;
    NM.TrainParam.LABEL.origlabel = NM.label;
    NM.TrainParam.SYNTH.flag    = 2;
elseif ~isfield(NM.TrainParam,'verbosity')
    NM.TrainParam.verbosity     = 1;
    if size(NM.label,2)>1
        NM.TrainParam.MULTILABEL.flag = true;
        NM.TrainParam.MULTILABEL.sel = 1:size(NM.label,2);
        NM.TrainParam.MULTILABEL.dim = size(NM.label,2);
    end
end

% Check whether the number of modalities specified in FUSION.M matches
% the number of modalities in NM (Users were copying their TrainParam
% structures from one NM workspace to the other and in cases where a data
% fusion was implemented in on NM workspace but was not possible in another
% workspace NM crashed.
if isfield(NM.TrainParam,'FUSION') && isfield(NM.TrainParam.FUSION,'M')
    nM1 = numel(NM.Y);
    nM2 = numel(NM.TrainParam.FUSION.M);
    if nM1 < nM2 || nM1 < max(NM.TrainParam.FUSION.M)
        NM.TrainParam.FUSION.M = 1;
        NM.TrainParam.FUSION.flag = 0;
    end
end
NM.TrainParam.ActiveModality = varind;

% for compatibility reasons
if ~isfield(NM.TrainParam,'STACKING') || ~NM.TrainParam.STACKING.flag
    NM.TrainParam.STACKING.flag=2;
end
if ~isfield(NM.TrainParam,'SYNTH') || ~NM.TrainParam.SYNTH.flag
    NM.TrainParam.SYNTH.flag=2;
end

% Adjust varind when in stacking mode and remove any SPATIAL operations
% from current modality
if NM.TrainParam.STACKING.flag == 1
    NM.TrainParam.FUSION.M = 1;
    NM.TrainParam.FUSION.flag=0;
    varind = 1;
    if isfield(NM.TrainParam.PREPROC{varind},'SPATIAL')
        NM.TrainParam.PREPROC{varind} = rmfield(NM.TrainParam.PREPROC{varind},'SPATIAL');
    end
    NM.TrainParam.ActiveModality = varind;
end

nan_in_label=false;         if sum(isnan(NM.label(:)))>0, nan_in_label=true; end
nY = numel(NM.Y);

if ~isfield(NM.TrainParam,'MLI')
    NM.TrainParam.MLI=[];
    for i=1:nY
        NM.TrainParam.MLI = nk_MLI_config(NM.TrainParam.MLI, i, 1);
    end
end

%% Create further default configurations
if ~isfield(NM.TrainParam,'PREPROC')
    % Create PREPROC structure
    for i=1:nY
        nan_in_pred = false;        if sum(isnan(NM.Y{i}(:)))>0, nan_in_pred=true; end
        NM.TrainParam.PREPROC{i}    = DefPREPROC(modeflag,nan_in_pred,nan_in_label);
        NM.TrainParam.VIS{i}        = nk_Vis_config([], NM.TrainParam.PREPROC, i, 1);
    end
else
    switch NM.TrainParam.FUSION.flag
        case 3 % If late fusion has been activated create STRAT structure
            for i=1:nY
                if isempty(NM.TrainParam.STRAT{i})
                    NM.TrainParam.STRAT{i}.SVM          = NM.TrainParam.SVM;
                    NM.TrainParam.STRAT{i}.GRD          = NM.TrainParam.GRD;
                    if numel(NM.TrainParam.PREPROC) == i
                        NM.TrainParam.STRAT{i}.PREPROC  = NM.TrainParam.PREPROC{i};
                        NM.TrainParam.STRAT{i}.VIS      = NM.TrainParam.VIS{i};
                    else
                        NM.TrainParam.STRAT{i}.PREPROC  = NM.TrainParam.PREPROC{1};
                        NM.TrainParam.STRAT{i}.VIS      = NM.TrainParam.VIS{1};
                    end
                    NM.TrainParam.STRAT{i}.RFE          = NM.TrainParam.RFE;
                    NM.TrainParam.STRAT{i}.MULTI        = NM.TrainParam.MULTI;
                    NM.TrainParam.STRAT{i}.MLI          = NM.TrainParam.MLI;
                    NM.TrainParam.STRAT{i}.MLI.Modality{1} = NM.TrainParam.MLI.Modality{i};
                end
            end
        otherwise
            nY = numel(NM.Y);
            nP = numel(NM.TrainParam.PREPROC);
            for i=1:nY
                nan_in_pred=false;          if sum(isnan(NM.Y{i}(:)))>0, nan_in_pred=true; end
                if ~isfield(NM.TrainParam,'PREPROC') || numel(NM.TrainParam.PREPROC) < nP
                    NM.TrainParam.PREPROC{i}    = DefPREPROC(modeflag, nan_in_pred, nan_in_label);
                end
                if ~isfield(NM.TrainParam,'VIS') || numel(NM.TrainParam.VIS) < nP
                    NM.TrainParam.VIS{i}        = nk_Vis_config([], NM.TrainParam.PREPROC{i}, i, 1);
                end
            end
    end
end

%% Check data entry status
if isfield(NM.TrainParam,'FUSION') && NM.TrainParam.FUSION.flag == 3
    STATUS = nk_CheckFieldStatus(NM,{'TrainParam','cv'},{'RAND', 'SYNTH', 'SAV', 'OOCV', 'META', 'STACKING', 'CALIB', 'LABEL'});
    STATUS = nk_CheckFieldStatus(NM.TrainParam.STRAT{varind},{'PREPROC','SVM','GRD','RFE','MULTI','VIS','MLI'}, [], [], STATUS);
else
    STATUS = nk_CheckFieldStatus(NM,{'TrainParam','cv'},{'STACKING','RAND', 'SYNTH','PREPROC','SVM','GRD','RFE','MULTI','VIS','SAV','OOCV','MLI', 'CALIB', 'LABEL'});
end
switch STATUS.PREPROC
    case '...'
        STATUS.FEATGEN = '...';
    otherwise
        STATUS.FEATGEN = '???';
end

switch NM.TrainParam.verbosity
    case 0
        verbostr = 'No output';
    case 1
        verbostr = 'Detailed output';
end

if ~exist('act','var') || isempty(act)

    %% Check whether there are analyses that have been completed and make stacking options available
    s = nk_GetNMStatus(NM);
    if ~isempty(s.completed_analyses) && sum(s.completed_analyses)>1 && sum(s.nmodal_analyses)>1
        menustr = [ menustr sprintf('Define meta-learning/stacking options [ %s ]|', STATUS.STACKING) ]; menuact = [ menuact 19 ];
    else
        NM.TrainParam.STACKING.flag = 2;
    end

    %% Check whether more than one variate are available and make data fusion options available
    if length(NM.Y)>1 && (NM.TrainParam.STACKING.flag == 2 || (~isfield(NM,'analysis') || (isfield(NM,'analysis') && numel(NM.analysis)<2)))
        % Make data fusion option available
        if isfield(NM.TrainParam,'FUSION')
            fusemode = NM.TrainParam.FUSION.flag;
            switch fusemode
                case 0
                    fusstr = ['[ Disabled. All operations will be applied to modality #' num2str(varind) ' ]'];
                case 1
                    fusstr  = '[ BEFORE preprocessing => ';
                case 2
                    fusstr  = '[ AFTER preprocessing => ';
                case 3
                    fusstr  = '[ Decision-based fusion (bagging) => ';
            end
            if fusemode
                fusesel = NM.TrainParam.FUSION.M;
                nF = numel(fusesel);
                if nF > 3
                    fusstr = [ fusstr sprintf('%g modalities ]',nF) ];
                else
                    fusstr = [ fusstr 'modalities:' ];
                    for i=1:nF, fusstr = sprintf('%s #%g +', fusstr, NM.TrainParam.FUSION.M(i)); end
                    fusstr = [fusstr(1:end-2) ' ]'];
                end
            end
        else
            fusstr = '[ undefined ]';
            fusemode = 0;
        end
        fusstr = ['Define data fusion options ' fusstr '|']; menustr = [menustr fusstr]; menuact = [ menuact 1 ];
        switch fusemode
            case {2,3}
                % Make active modality selection option available
                descstr = [ ' (' NM.datadescriptor{varind}.desc ')'];
                varstr = ['Set active modality for configuration [ #' num2str(varind) descstr ' ]|']; menuact = [ menuact 2 ];
                if fusemode == 3, NM.TrainParam.META.flag = 1; end
                resetstr = '';
                if fusemode == 3
                    resetstr = 'Reset all modality configurations to current modality setup|';  menuact = [ menuact 5 ];
                end
                menustr = [menustr varstr resetstr];
        end
    else
        fusemode=0;
    end

    % Check for regression / classification experiment and number of groups
    multistr = ''; multiflag = false;
    if isfield(NM.TrainParam, 'LABEL') && NM.TrainParam.LABEL.flag
        modeflag = NM.TrainParam.LABEL.newmode;
    else
        modeflag = NM.modeflag;
        NM.TrainParam.LABEL.flag = 0; 
    end

    if strcmp(modeflag,'classification')

        classtr = ['Classification algorithm [ ' STATUS.SVM ' ]|'];
        
        if isfield(NM.TrainParam, 'LABEL') && NM.TrainParam.LABEL.flag
            label_temp = NM.TrainParam.LABEL.newlabel;
        else
            label_temp = NM.label;
        end
        if numel(unique(label_temp(~isnan(label_temp)))) > 2
            multistr = ['Multi-class settings [ ' STATUS.MULTI ' ]|'];
            multiflag = true;
        else
            NM.TrainParam.MULTI.flag=0;
        end
    else
        classtr = ['Prediction algorithm [ ' STATUS.SVM ' ]|'];
        multistr = ''; NM.TrainParam.MULTI.flag=0;
    end

    if isfield(NM,'OOCV') && ~isempty(NM.OOCV)
        oocvstr = ['OOCV prediction parameters [ ' STATUS.OOCV ' ]|'];
        oocvflag = true;
    else
        oocvstr='';
        oocvflag = false;
    end
    if size(NM.label,2)>1
        if isfield(NM,'TrainParam') && isfield(NM.TrainParam,'MULTILABEL')
            nL = nk_GetLabelDim(NM.TrainParam.MULTILABEL);
            if nL>1
                multlabelselstr = sprintf('%g labels selected',nL);
            else
                multlabelselstr = sprintf('Label: %s',NM.labelnames{NM.TrainParam.MULTILABEL.sel});
            end
        else
            multlabelselstr = sprintf('%g labels selected',size(NM.label,2));
        end
        menustr = [ menustr 'Label selection in multi-label mode [ ' multlabelselstr ' ]|']; menuact = [ menuact 21];
    end

    menustr = [menustr 'Cross-validation settings [ ' STATUS.cv  ' ]|']; menuact = [menuact 3];

    SVM = [];
    if isfield(NM.TrainParam,'FUSION') && NM.TrainParam.FUSION.flag == 3
        flFUS = true;
        flSVM = isfield(NM.TrainParam.STRAT{varind},'SVM');
        flGRD = isfield(NM.TrainParam.STRAT{varind},'GRD');
        flPREPROC = isfield(NM.TrainParam.STRAT{varind},'PREPROC');
        if flSVM, SVM = NM.TrainParam.STRAT{varind}.SVM; end
    else
        flFUS = false;
        flSVM = isfield(NM.TrainParam,'SVM');
        flGRD = isfield(NM.TrainParam,'GRD');
        flPREPROC = isfield(NM.TrainParam,'PREPROC');
        if flSVM, SVM = NM.TrainParam.SVM; end
    end

    flx = flSVM && flGRD && flPREPROC;

    menustr = [ menustr 'Use different label [ ' STATUS.LABEL ' ]|']; menuact = [menuact 99] ;
    
    menustr = [ menustr sprintf('Synthethic data generation [ %s ]|', STATUS.SYNTH) ]; menuact = [ menuact 4 ];

    menustr = [ menustr 'Preprocessing pipeline [ ' STATUS.PREPROC ' ]|' classtr ]; menuact = [ menuact 6:7 ];

    if flx

        if (any(strcmp(SVM.prog,{'MikRVM','MKLRVM','MVTRVR','MVTRVM','GLMFIT'})) && any(strcmp(SVM.kernel.kernstr,{' -t 0','lin','linear'}))) || ...
                ( strcmp(SVM.prog,'matLRN') && ( ~isfield(NM.TrainParam.GRD.matLearn,'Params') || isempty(NM.TrainParam.GRD.matLearn.Params)) )
            if flFUS
                NM.TrainParam.STRAT{varind}.GRD.PX = []; NM.TrainParam.STRAT{varind}.GRD.n_params = 0;
            else
                NM.TrainParam.GRD.PX = []; NM.TrainParam.GRD.n_params = 0;
            end
        else
            menustr = [ menustr 'Learning algorithm parameters [ ' STATUS.GRD ' ]|' ]; menuact = [ menuact 8 ];
        end
        % For the predictive sequence optimizer we don't need any ensemble
        % learning optimization.
        if ~strcmp(SVM.prog,'SEQOPT')
            menustr = [ menustr ...
                'Ensemble generation strategies [ ' STATUS.RFE ' ]|'];
            menuact = [menuact 9];
        end
        if multiflag
            menustr = [ menustr multistr];
            menuact = [ menuact 10 ];
        end
    end

    menustr = [ menustr ...
        'Visualization options [ ' STATUS.VIS ' ]|' ...
        'Prediction interpretation options [ ' STATUS.MLI ' ]|' ...
        'Model saving options [ ' STATUS.SAV ' ]|' ...
        oocvstr ...
        'Define verbosity level [ ' verbostr ' ]|' ...
        'Inspect workspace|' ...
        'Save parameter template|' ...
        'Load training template'];


    menuact = [ menuact 11 12 13 ];
    if oocvflag, menuact = [ menuact 14 ]; end
    menuact = [ menuact 17 20 15 16 ];

    nk_PrintLogo
    mestr = 'Define parameter template'; navistr = sprintf('%s\n\t>>> %s',parentstr, mestr); fprintf('\nYou are here: %s >>> ',parentstr);
    if fusemode > 1, fprintf('\n==> CONFIGURATION OF MODALITY #%g: %s', varind, descstr ); end
    act = nk_input(mestr, 0, 'mq', menustr, menuact);
end

switch act

    case 0
        return
    case 1
        if isfield(NM.TrainParam,'FUSION')
            fusedef = NM.TrainParam.FUSION.flag ;
        else
            fusedef = 1;
        end
        NM.TrainParam.FUSION.flag = nk_input('Select multimodal fusion strategy',0,'m', ...
            ['No fusion|' ...
            'Early fusion -> Modality concatenation BEFORE feature preprocessing|' ...
            'Intermediate fusion -> Modality concatenation AFTER preprocessing|' ...
            'Late fusion -> Decision-based data fusion (bagging)'],0:3, fusedef);
        if NM.TrainParam.FUSION.flag
            NM.TrainParam.FUSION.M = nk_SelectVariateIndex(NM, 0, 1, 1 );
            if numel(NM.TrainParam.FUSION.M) == 1
                NM.TrainParam.FUSION.flag = false;
            end
            if isempty(find(NM.TrainParam.FUSION.M == varind, 1))
                varind = NM.TrainParam.FUSION.M(1);
            end
            if NM.TrainParam.FUSION.flag == 3
                if ~isfield(NM.TrainParam,'STRAT')
                    NM.TrainParam.STRAT = cell(numel(NM.Y),1);
                end
            end
        else
            % Only one modality has to be selected for analysis
            NM.TrainParam.FUSION.M = nk_SelectVariateIndex(NM, 1, 1, 1 );
            varind = NM.TrainParam.FUSION.M;
        end
        NM.TrainParam.ActiveModality = varind;

    case 2
        if isfield(NM.TrainParam,'FUSION') && NM.TrainParam.FUSION.flag 
            M = NM.TrainParam.FUSION.M ;
        else
            M = [];
        end
        varind = nk_SelectVariateIndex(NM, 1, 1, 1, M);
        NM.TrainParam.ActiveModality = varind;

    case 3
        act = 1; while act>0, act = nk_CVpartition_config; end
    
    case 4
        if ~isfield(NM.TrainParam,'SYNTH'), NM.TrainParam.SYNTH.flag = 2; end
        mess=[];act = 1; while act>0, [NM.TrainParam.SYNTH, act, mess] = nk_Synth_config(NM.TrainParam.SYNTH, mess, navistr); end

    case 5
        fl = questdlg('Are you sure you want to overwrite all modality configuration with the current modality setup?',mestr,'Yes','No','No');
        if strcmp(fl,'Yes')
            STRAT = NM.TrainParam.STRAT(varind);
            NM.TrainParam.STRAT(NM.TrainParam.FUSION.M) = STRAT;
        end

    % PREPROCESSING =================================================================================================================================================
    case 6
        if ~isfield(NM,'TrainParam') || ...
                ~isfield(NM.TrainParam,'PREPROC') || ...
                varind > numel(NM.TrainParam.PREPROC)
            NM.TrainParam.PREPROC{varind}.FEATSEL.active=0;
            NM.TrainParam.PREPROC{varind}.BINMOD=1;
        end
        if isfield(NM.TrainParam,'FUSION') && NM.TrainParam.FUSION.flag == 3
            if ~isfield(NM.TrainParam.STRAT{varind},'PREPROC'), NM.TrainParam.STRAT{varind}.PREPROC = []; end
            act = 1; stepind = 1; while act>0, [NM.TrainParam.STRAT{varind}.PREPROC, act, stepind] = nk_Preproc_config(NM.TrainParam.STRAT{varind}.PREPROC, varind, navistr, stepind); end
        else
            if isfield(NM.TrainParam,'FUSION') && NM.TrainParam.FUSION.flag == 1
                varind = NM.TrainParam.FUSION.M(1); %Always use the first modality to setup the parameters for the concatenated feature space
            end
            if ~isfield(NM.TrainParam,'PREPROC'), NM.TrainParam.PREPROC{varind} = []; end
            act = 1; stepind = 1; while act>0, [NM.TrainParam.PREPROC{varind}, act, stepind] = nk_Preproc_config(NM.TrainParam.PREPROC{varind}, varind, navistr, stepind); end
        end

        % ML ALGORITHM SELECTION =========================================================================================================================================
    case 7
        if isfield(NM.TrainParam,'FUSION') && NM.TrainParam.FUSION.flag == 3
            if ~isfield(NM.TrainParam.STRAT{varind},'SVM'), NM.TrainParam.STRAT{varind}.SVM = []; end
            act = 1; while act>0, ...
                    [NM.TrainParam.STRAT{varind}.SVM, act] = nk_Model_config(NM.TrainParam.STRAT{varind}.SVM, NM.TrainParam.STRAT{varind}, navistr, varind); end
        else
            if ~isfield(NM.TrainParam,'SVM'), NM.TrainParam.SVM = []; end
            act = 1; while act>0, [NM.TrainParam.SVM, act ] = nk_Model_config(NM.TrainParam.SVM, NM.TrainParam, navistr, varind); end
        end

        % ML OPTIMIZATION STRATEGIES =====================================================================================================================================
    case 8
        if isfield(NM.TrainParam,'FUSION') && NM.TrainParam.FUSION.flag == 3
            if ~isfield(NM.TrainParam.STRAT{varind},'GRD'), NM.TrainParam.STRAT{varind} = nk_Grid_config(NM.TrainParam.STRAT{varind}, NM.TrainParam.STRAT{varind}.SVM, [], true); end
            act = 1; while act>0, ...
                    [ NM.TrainParam.STRAT{varind}, act ] = nk_Grid_config(NM.TrainParam.STRAT{varind}, NM.TrainParam.STRAT{varind}.SVM, [], [], navistr); end
        else
            if ~isfield(NM.TrainParam,'GRD'), NM.TrainParam = nk_Grid_config([], NM.TrainParam.SVM, varind, true); end
            act = 1; while act>0, [ NM.TrainParam, act ] = nk_Grid_config(NM.TrainParam, NM.TrainParam.SVM, varind, [], navistr); end
        end
        % FEATURE SELECTION ==============================================================================================x=================================================
    case 9
        if isfield(NM.TrainParam,'FUSION') && NM.TrainParam.FUSION.flag == 3
            if ~isfield(NM.TrainParam.STRAT{varind},'RFE')
                [~, NM.TrainParam.STRAT{varind}.RFE ] = ...
                    nk_RFE_config([], NM.TrainParam.STRAT{varind}, ...
                    NM.TrainParam.STRAT{varind}.SVM, ...
                    modeflag, ...
                    NM.TrainParam.STRAT{varind}.MULTI, ...
                    NM.TrainParam.STRAT{varind}.GRD, 1);
            end
            act = 1;
            while act>0
                [ act, NM.TrainParam.STRAT{varind}.RFE ] = ...
                    nk_RFE_config(act, NM.TrainParam.STRAT{varind}, ...
                    NM.TrainParam.STRAT{varind}.SVM, ...
                    modeflag, ...
                    NM.TrainParam.STRAT{varind}.MULTI, ...
                    NM.TrainParam.STRAT{varind}.GRD, [], navistr);
            end
        else
            if ~isfield(NM.TrainParam,'RFE'), [~, NM.TrainParam.STRAT{varind}.RFE ] = ...
                    nk_RFE_config([], NM.TrainParam, NM.TrainParam.SVM, modeflag, NM.TrainParam.MULTI, NM.TrainParam.GRD, 1); end
            act = 1; while act>0, [ act, NM.TrainParam.RFE ] = nk_RFE_config(act, NM.TrainParam, NM.TrainParam.SVM, modeflag, NM.TrainParam.MULTI, NM.TrainParam.GRD, [], navistr); end
        end

        % MULTI-GROUP SETTINGS =============================================================================================================================================
    case 10
        if isfield(NM.TrainParam,'FUSION') && NM.TrainParam.FUSION.flag == 3
            if ~isfield(NM.TrainParam.STRAT{varind},'MULTI'), NM.TrainParam.STRAT{varind}.MULTI = nk_Multi_config([], true); end
            act = 1; while act>0, [ NM.TrainParam.STRAT{varind}.MULTI, act ] = nk_Multi_config(NM.TrainParam.STRAT{varind}.MULTI, [], navistr); end
        else
            if ~isfield(NM.TrainParam,'MULTI'), NM.TrainParam.MULTI = nk_Multi_config([], true); end
            act = 1; while act>0, [ NM.TrainParam.MULTI, act] = nk_Multi_config(NM.TrainParam.MULTI,[], navistr); end
        end

        % VISUALIZATION ====================================================================================================================================================
    case 11
        if isfield(NM.TrainParam,'FUSION') && NM.TrainParam.FUSION.flag == 3
            if ~isfield(NM.TrainParam,'VIS'), NM.TrainParam.STRAT{varind}.VIS = nk_Vis_config(NM.TrainParam.STRAT{varind}.VIS, NM.TrainParam.STRAT{varind}.PREPROC, 1, 1, navistr); end
            act = 1; while act>0, [ NM.TrainParam.STRAT{varind}.VIS, act] = nk_Vis_config(NM.TrainParam.STRAT{varind}.VIS, NM.TrainParam.STRAT{varind}.PREPROC, 1, [], navistr); end
        else
            if ~isfield(NM.TrainParam,'VIS'), NM.TrainParam.VIS{varind} = nk_Vis_config(NM.TrainParam.VIS{varind}, NM.TrainParam.PREPROC{varind}, varind, 1, navistr); end
            act = 1; while act>0, [ NM.TrainParam.VIS{varind}, act ] = nk_Vis_config(NM.TrainParam.VIS{varind}, NM.TrainParam.PREPROC{varind}, varind , [], navistr); end
        end

        % ML INTEPRETATION STRATEGIES =======================================================================================================================================
    case 12
        if isfield(NM.TrainParam,'FUSION') && NM.TrainParam.FUSION.flag == 3
            if ~isfield(NM.TrainParam,'MLI'), NM.TrainParam.STRAT{varind}.MLI = nk_MLI_config(NM.TrainParam.STRAT{varind}.MLI, 1, 1, navistr); end
            act = 1; while act>0, [ NM.TrainParam.STRAT{varind}.MLI, act] = nk_MLI_config(NM.TrainParam.STRAT{varind}.MLI, 1, [], navistr); end
        else
            if ~isfield(NM.TrainParam,'MLI'), NM.TrainParam.MLI = nk_MLI_config([], varind, 1, navistr); end
            act = 1; while act>0, [ NM.TrainParam.MLI, act ] = nk_MLI_config(NM.TrainParam.MLI, varind, [], navistr); end
        end

        % SAVING OPTIONS ====================================================================================================================================================
    case 13
        act = 1; while act>0, [ NM, act ] = nk_SavingOptions_config(NM, 0, navistr); end

        % OOCV OPTIONS ======================================================================================================================================================
    case 14
        NM.TrainParam = nk_OOCV_config(NM.TrainParam);

        % EXPORT TRAINING PARAM =============================================================================================================================================
    case 15
        matname = nk_input('Filename (prefix is TRAIN)',0,'s');
        matname = ['TRAIN_' matname];
        if isfield(NM,'TrainParam')
            TrainParam = NM.TrainParam;
            save(matname,'TrainParam');
        end
        if isfield(NM,'cv')
            cv = NM.cv;
            save(matname,'cv','-append');
        end

        % IMPORT TRAINING PARAM =============================================================================================================================================
    case 16
        fl = nk_input('Overwrite current settings?',0,'yes|no',[1,0],1);
        if fl
            menuvec = []; menustr =[];
            if isfield(NM,'analysis')
                menustr = 'Take parameters from analysis in current NM structure|';
                menuvec = 1;
            end
            menustr = [menustr 'Load parameters from parameter file|'];
            menuvec(end+1) = 2;
            menustr = [menustr 'Load parameters from NM structure file'];
            menuvec(end+1) = 3;
            act = nk_input('Source of parameters to be loaded', 0, 'm', menustr, menuvec);

            switch act
                case {2,3}
                    mess = [];
                    switch act
                        case 2
                            matname = nk_FileSelector(1,'matrix','select TRAIN mat','TRAIN.*\mat');
                            if exist(matname,'file'),load(matname); end
                        case 3
                            matname = nk_FileSelector(1,'matrix','Select NM structure file','.*\mat');
                            if exist(matname,'file')
                                [~, matfile] = fileparts(matname);
                                fprintf('\nLoading %s as temporary structure',matfile)
                                load(matname,'NM','TrainParam');
                                load(matname,'NM','cv');
                            end
                    end
                    if exist('TrainParam','var')
                        if isfield(NM,'TrainParam') && isfield(NM.TrainParam,'RAND')
                            RAND = NM.TrainParam.RAND;
                        end
                        NM.TrainParam = TrainParam; mess{1} = 'TrainParam';
                        if exist('RAND','var'), NM.TrainParam.RAND = RAND; end
                    end
                    if exist('cv','var') && ~isfield(NM,'cv')
                        NM.cv = cv; mess{1+end} = 'CV';
                    else
                        mess{1+end} = 'CV already existed, not overwritten!';
                    end
                case 1
                    %[~,analind] = nk_SelectAnalysis(NM); mess=[];
                    t_act = 1; brief = 1; while t_act>0, [t_act, analind, NM, ~, brief ]= nk_SelectAnalysis(NM, 0, 'MAIN >>> Select Analysis', 1, 1,0,[],brief); end
                    if isfield(NM.analysis{analind}.params,'TrainParam')
                        mess{1} = 'TrainParam';
                        NM.TrainParam = NM.analysis{analind}.params.TrainParam;
                    end
                    if isfield(NM.analysis{analind}.params,'cv')
                        NM.cv = NM.analysis{analind}.params.cv;
                        mess{1+end} = 'CV';
                    end
                    NM.cv = NM.analysis{analind}.params.cv;
            end
            if ~isempty(mess)>0, msgbox(mess,'Loaded parameters:','none'); end
        end

        % DEFINE VERBOSITY LEVEL ============================================================================================================================================
    case 17
        NM.TrainParam.verbosity = ~NM.TrainParam.verbosity;

    case 18

        act = 1; while act>0, [NM.TrainParam.META, act] = nk_Ensemble_config(NM); end
        
        % META-LEARNING (STACKING) =======================================================================================================================================
    case 19
        if ~isfield(NM.TrainParam,'STACKING'), NM.TrainParam.STACKING.flag = 2; end
        mess=[];act = 1; while act>0, [NM.TrainParam.STACKING, act, mess] = nk_Stacking_config(NM.TrainParam.STACKING, s, mess, navistr); end

    case 20
        nk_PrintWs(NM, NM.TrainParam)

    case 21
        nk_PrintLogo
        fprintf('\n*************************************')
        fprintf('\n***     MULTI-LABEL SELECTION     ***')
        fprintf('\n*************************************')
        fprintf('\n')
        
        for i=1:size(NM.label,2)
            fprintf('\n- %g -  %s',i, NM.labelnames{i})
        end
        if ~isfield(NM.TrainParam,'MULTILABEL')
            NM.TrainParam.MULTILABEL.sel = 1:size(NM.label,2);
        end
        NM.TrainParam.MULTILABEL.sel = nk_input('Select labels for processing',0,'i',NM.TrainParam.MULTILABEL.sel);
        NM.TrainParam.MULTILABEL.dim = numel(NM.TrainParam.MULTILABEL.sel);
        if isfield(NM.TrainParam.RAND,'Eq') && NM.TrainParam.RAND.Eq.enabled && strcmp(NM.TrainParam.RAND.Eq.CovarName,'NM.label')
            NM.TrainParam.RAND.Eq.Covar = nm_nanmean(NM.label(:,NM.TrainParam.MULTILABEL.sel),2);
            nk_PrintLogo
            fprintf('\nFound active histogram equalization based on NM.label!')
            reinitcv = nk_input('Shall I reinitialize the CV structure ?', 0, 'yes|no', [1,0], 1);
            if reinitcv
                nk_CVpartition_config(0,6);
            end
        end

    case 99
        nk_PrintLogo
        fprintf('\n*************************************')
        fprintf('\n****  DEFINE ALTERNATIVE LABEL  *****')
        fprintf('\n*************************************')
        fprintf('\n')
        if isfield(NM.TrainParam, 'LABEL')
            LABEL = NM.TrainParam.LABEL;
        else
            LABEL = [];
        end
        if ~isfield(LABEL, 'OrigTrainParam')
            LABEL.OrigTrainParam = NM.TrainParam; 
        end
        if ~isfield(LABEL, 'origlablabel')
            LABEL.origlabel = NM.label;
        end
        
        while act>0  
            [LABEL, act] = cv_Label_config(LABEL);
        end
        if LABEL.flag && ~isempty(LABEL.newmode) && ~strcmp(LABEL.newmode, modeflag)
            origmodefl                  = NM.modeflag;
            % check whether a new mode was entered
            if isempty(LABEL.newmode)
                LABEL.newmode           = origmodefl;
            end
            NM.modeflag                 = LABEL.newmode;
            
            % Create default NM parameters space
            nk_CVpartition_config(true);
            NM.TrainParam.STACKING.flag = 2;
            NM.TrainParam.FUSION.flag   = 0;
            NM.TrainParam.FUSION.M      = 1;
            [~,NM.TrainParam.SVM]           = nk_LIBSVM_config(NM,[],1);
            NM.TrainParam.SVM.prog      = 'LIBSVM';
            NM.TrainParam.SVM           = nk_Kernel_config(NM.TrainParam.SVM,1);
            NM.TrainParam.SVM.GridParam = 1;
            if strcmp(NM.modeflag, 'regression'), NM.TrainParam.SVM.GridParam = 18; end
            NM.TrainParam.MULTI.flag    = 0;
            NM.TrainParam               = nk_Grid_config(NM.TrainParam, NM.TrainParam.SVM, varind, true);
            [~,NM.TrainParam.RFE]       = nk_RFE_config([], NM.TrainParam, NM.TrainParam.SVM, modeflag, NM.TrainParam.MULTI, NM.TrainParam.GRD, 1);
            NM.TrainParam.verbosity     = 1;

            NM.TrainParam               = rmfield(NM.TrainParam,'PREPROC');

            NM.TrainParam.LABEL         = LABEL;
            NM.modeflag                 = origmodefl;

        elseif LABEL.flag % but same learning framework
            NM.TrainParam.LABEL         = LABEL;
        elseif ~LABEL.flag && strcmp(LABEL.newmode, modeflag) % if switched from alternative label to no alt. label but no learning mode switch
            NM.TrainParam = rmfield(NM.TrainParam, 'LABEL'); 
        elseif ~LABEL.flag % if either nothing has been changed in submenu or switch from alternative label to no alt. label and learning mode switch 
            NM.TrainParam = LABEL.OrigTrainParam;
        end
      


 %% read in calibration data
    case 1000

        NM = nk_DefineOOCVData_config(NM, 2, 'calib');
        NM = nk_SelectOOCVdata(NM, 2, 0);
        CALIBAVAIL = 1;
        CY = NM.C{1,1}.Y;
        Cfile_path = sprintf('%s/CY.mat', pwd);
        NM.C{1,1}.calibflag = 1;
        save(Cfile_path, 'CY', '-v7.3');
        NM.C{1,1}.Y = Cfile_path;
end
act = 1;

function PREPROC = DefPREPROC(modeflag, nan_in_pred, nan_in_label)

if ~exist('nan_in_pred','var') || isempty(nan_in_pred), nan_in_pred = false; end
if ~exist('nan_in_label','var') || isempty(nan_in_label), nan_in_label = false; end

PREPROC.BINMOD = 1; PREPROC.FEATSEL.active = 0;

PREPROC.ACTPARAM{1}.SCALE   = nk_Scale_config([],[],1);
PREPROC.ACTPARAM{1}.cmd     = 'scale';

PREPROC.ACTPARAM{2}.PRUNE   = nk_Prune_config([],[],[],1);
PREPROC.ACTPARAM{2}.cmd     = 'elimzero';

if nan_in_pred % Adjust defaults to NaN in the predictor data
    PREPROC.ACTPARAM{1}.SCALE.zerooutflag   = 1;
    PREPROC.ACTPARAM{2}.PRUNE.nan           = 2;
    PREPROC.ACTPARAM{2}.PRUNE.inf           = 2;
    PREPROC.ACTPARAM{3}.IMPUTE              = nk_Impute_config([],[],[],[],1);
    PREPROC.ACTPARAM{3}.cmd                 = 'impute';
end

if nan_in_label % Adjust defaults to NaN in the labels data
    PREPROC.LABELMOD.LABELIMPUTE            = nk_Impute_config([], [], [], [], 1);
    PREPROC.LABELMOD.cmd                    = 'labelimpute';
end

if strcmp(modeflag,'regression'), PREPROC.LABELMOD.TARGETSCALE = true; end
