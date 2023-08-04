function [ act, NM, inp ] = cv_ExportModelsPrep(NM, act, inp, parentstr)
global CV

% Detect completed analyses
as = nk_GetAnalysisStatus(NM); complvec = find(as.completed_analyses);

% Initialize runtime parameters
if ~exist('inp','var') || isempty(inp)
    inp = struct( 'analind', complvec(1), ...   % Index to analysis
        'oocvind', 1, ...           % Index to OOCV data container                    
        'lfl', 1, ...               % 1 = compute from scratch |
        ...                         % 2 = use existing (allowing the user to specify OOCVdatamats)
        'ovrwrt', 2, ...            % if lfl == 1 ==> 1 = overwrite existing OOCVdatamats
        ...                         % 2 = do not overwrite (use existing OOCVdatamats automatically)
        'saveparam', 1, ...         % if loadparam == 2=> 1 = save OOCV processing parameters (preprocessing / models)
        ...                         % 2 = do not save parameters to disk
        'saveCV1', 2, ...           % if loadparam == 2 && saveparam ==1 => 1 = save large OOCV processing containers at the CV1 level
        ...                         % 2 = operate at CV2 level
        'loadparam', 2, ...         % 1 = load existing optpreproc and/or optmodel parameters from disk
        ...                         % 2 = recompute parameters
        'HideGridAct', false, ...
        'batchflag', 0);            % 1 = Run in batchmode (without graphics outputs)
    % 0 = run in interactive mode
end

na_str = '?'; 
OverWriteStr = []; GridSelectStr = []; LoadModelsStr = []; LoadParamsStr = []; LoadStr = []; SaveStr = []; SaveCV1Str = [];
OverWriteAct = []; GridSelectAct = []; LoadModelsAct = []; LoadParamsAct = []; LoadAct = []; SaveAct = []; SaveCV1Act = [];

%% Configure menu
% Select analysis
if numel(NM.analysis)>1
    if numel(inp.analind)<2
        AnalSelStr = sprintf('Analysis %g', inp.analind);
    else
        if ~inp.HideGridAct, cvequalstr = 'same-size CV structures'; else, cvequalstr = 'different CV structures'; end
        AnalSelStr = sprintf('%g Analyses: %s [ %s ]',numel(inp.analind), strjoin(cellstr(num2str(inp.analind'))',', '), cvequalstr);
    end
    AnalSelectStr = sprintf('Choose analysis to work on [ %s ]|', AnalSelStr);                                              AnalSelectAct = 1;
else
    AnalSelectStr = ''; AnalSelectAct = [];
end

analysis      = NM.analysis{inp.analind(1)};

if ~isempty(analysis)

    % Initialize global parameters for the selected analysis
    nk_SetupGlobalVariables(analysis.params, 'setup_main', 0);

    % Retrieve CV2 partitions to operate on
    if ~isfield(inp,'GridAct'), inp.GridAct = analysis.GDdims{1}.GridAct; end
    if ~inp.HideGridAct
        GridSelectStr = sprintf('Select CV2 partitions to operate on [ %g selected ]|',  sum(inp.GridAct(:)));              GridSelectAct = 5;
    else
        GridSelectStr =''; GridSelectAct=[];
    end
    if inp.loadparam == 2 && inp.lfl == 1
        SAVE_opts       = {'yes', 'no'};
        SaveStr = sprintf('Save pre-processing params and models to disk [ %s ]|', SAVE_opts{inp.saveparam});               SaveAct = 6;
        if inp.saveparam == 1
            SaveCV1Str = sprintf('Save pre-processing params at CV1 level [ %s ]|', SAVE_opts{inp.saveCV1});                SaveCV1Act = 12;
        end
    end
end

%% Build interactive menu
menustr = [ AnalSelectStr ...
    OverWriteStr ...
    GridSelectStr ...
    ];

menuact = [ AnalSelectAct ...
    OverWriteAct ...
    GridSelectAct ... 
    ];

disallow = false;

%% Check whether all parameters are available
if (~sum(inp.GridAct(:)) && ~inp.HideGridAct) || isempty(inp.analind), disallow = true; end

if inp.loadparam == 1
    if ~isfield(inp,'optpreprocmat') || isempty(inp.optpreprocmat), disallow = true; end
    if ~isfield(inp,'optmodelmat') || isempty(inp.optmodelmat), disallow = true; end
end

if ~disallow, menustr = [menustr '|PROCEED >>>']; menuact = [menuact 10]; end

%% Display menu and act on user selections
nk_PrintLogo


mestr = 'Export model parameters'; navistr = [parentstr ' >>> ' mestr]; fprintf('\nYou are here: %s >>>',parentstr);
if ~inp.batchflag && act<13, act = nk_input(mestr, 0, 'mq', menustr, menuact); end

    switch act

        case 0
            return
            % Select analysis
        case 1
            showmodalvec = []; analind = inp.analind;
            if length(NM.analysis)>1, t_act = 1; brief = 1;
                while t_act>0
                    [t_act, analind, ~, showmodalvec , brief] = nk_SelectAnalysis(NM, 0, navistr, analind, [], 1, showmodalvec, brief);
                end
                if ~isempty(analind), inp.analind = complvec(analind) ; end
                nA = numel(inp.analind);
                if nA>1
                    AS = nk_GetAnalysisStatus(NM, inp.analind);
                    if ~AS.betweenfoldpermequal_cv
                        inp.HideGridAct = true;
                    else
                        inp.GridAct = NM.analysis{inp.analind(1)}.GDdims{1}.GridAct;
                        inp.HideGridAct = false;
                    end
                else
                    inp.HideGridAct = false;
                    inp.GridAct = NM.analysis{inp.analind}.GDdims{1}.GridAct;
                end
            end
        case 4
            switch inp.lfl
                case 1
                    if inp.ovrwrt == 1, inp.ovrwrt = 2; elseif inp.ovrwrt  == 2, inp.ovrwrt = 1; end
                case 2
                    tdir = create_defpath(NM.analysis{inp.analind}, inp.oocvind); %%TO DO EXPORTMODELS NECESSARY???
            end
        case 5
            [operms,ofolds] = size(CV.TrainInd);
            tact = 1; while tact > 0 && tact < 10, [ tact, inp.GridAct ] = nk_CVGridSelector(operms, ofolds, inp.GridAct, 0); end
        case 6
            if inp.saveparam == 1, inp.saveparam = 2; elseif inp.saveparam == 2,  inp.saveparam = 1; end
        case 7
            if inp.loadparam == 1, inp.loadparam = 2; elseif inp.loadparam == 2,  inp.loadparam = 1; end
        case 8
            tdir = create_defpath(NM.analysis{inp.analind}, inp.oocvind(1)); %%TO DO EXPORTMODELS NECESSARY???
            optpreprocmat = nk_GenDataMaster(NM.id, 'OptPreprocParam', CV, [], tdir);
            if ~isempty(optpreprocmat), inp.optpreprocmat = optpreprocmat; end
        case 9
            tdir = create_defpath(NM.analysis{inp.analind}, inp.oocvind(1)); %%TO DO EXPORTMODELS NECESSARY???
            optmodelmat = nk_GenDataMaster(NM.id, 'OptModel', CV, [], tdir);
            if ~isempty(optmodelmat), inp.optmodelmat = optmodelmat; end
        case {10,11}
            for j=1:numel(inp.oocvind) %%TO DO EXPORTMODELS NECESSARY???
                inp.oocvname = sprintf('OOCV_%g',inp.oocvind(j));
            end
            nA = 1; if numel(inp.analind)>1, nA = numel(inp.analind); end
            for i=1:nA
                NM.runtime.curanal = inp.analind(i);
                inp = nk_GetAnalModalInfo_config(NM, inp);
                if inp.HideGridAct, [ ix, jx ] = size(NM.analysis{inp.analind(i)}.params.cv.TrainInd); inp.GridAct = true(ix,jx); end
                inp.analysis_id = NM.analysis{inp.analind(i)}.id;
                inp.saveoptdir = [ NM.analysis{inp.analind(i)}.rootdir filesep 'opt' ];
                ExportModelsPrep(NM, inp, NM.analysis{inp.analind(i)});
                nk_SetupGlobalVariables(NM.analysis{inp.analind(i)}.params, 'clear', 0);
            end
            NM = rmfield(NM,'runtime');
        case 12
            if inp.saveCV1 == 1, inp.saveCV1 = 2; elseif inp.saveCV1 == 2,  inp.saveCV1 = 1; end
    end

    act = 1; 

end

function tdir = create_defpath(analysis, oocvind)

rootdir = analysis.GDdims{1}.RootPath;
if isfield(analysis,'OOCV') && numel(analysis.OOCV) >= oocvind && isfield(analysis.OOCV{oocvind},'RootPath')
    if iscell(analysis.OOCV{oocvind}.RootPath)
        tdir = analysis.OOCV{oocvind}.RootPath{1};
    else
        tdir = analysis.OOCV{oocvind}.RootPath;
    end
else
    oocvdir = sprintf('OOCV_%g', oocvind);
    tdir = fullfile(rootdir, oocvdir);
end
end
%
% =========================================================================
function ExportModelsPrep(dat, inp1, analysis)
global SAV MODEFL CV FUSION MULTILABEL
% tOOCV = OOCV;
if inp1.saveparam   == 2, inp1.saveparam    = 0; end
if inp1.ovrwrt      == 2, inp1.ovrwrt       = 0; end
if inp1.lfl         == 1, inp1.analmode     = 0; else, inp1.analmode = 1; end

F = 1; nF = 1;
if ~isempty(FUSION)
    F = analysis.params.TrainParam.FUSION.M;
    nF = numel(F); if FUSION.flag < 3, nF = 1; end
    inp1.nF = nF;
end

if strcmp(MODEFL,'classification')
    inp1.nclass = length(CV.class{1,1});
else
    inp1.nclass = 1;
end

inp1.id              = dat.id;
stranalysis          = SAV.matname;
inp1.ngroups         = numel(unique(dat.label(~isnan(dat.label))));

if isfield(inp1,'targdir') %%TO DO EXPORTMODELS NECESSARY???
    inp1.rootdir = fullfile(inp1.targdir, inp1.oocvname);
elseif isfield(analysis,'rootdir') && exist(analysis.rootdir,'dir')
    inp1.rootdir = fullfile(analysis.rootdir,analysis.params.TrainParam.SVM.prog, inp1.oocvname);
else
    inp1.rootdir = fullfile(pwd,analysis.params.TrainParam.SVM.prog, inp1.oocvname);
end

if ~exist(inp1.rootdir,'dir'), mkdir(inp1.rootdir); end
nl = nk_GetLabelDim(MULTILABEL);

% Loop through modalities
for i = 1:inp1.nF

    % **************************** ANALYSIS SETUP *****************************
    inp2 = nk_DefineFusionModeParams(dat, analysis, F, nF, i);
    inp = catstruct(inp1,inp2);
    inp.loadGD = true;

    for j = 1:nl

        inp.multlabelstr = '';  if MULTILABEL.flag, inp.multlabelstr = sprintf('_t%g',j); end

        if MULTILABEL.flag && MULTILABEL.dim>1
            fprintf('\n\n');fprintf('====== Working on label #%g ====== ',j);
            inp.curlabel = j;
        else
            inp.curlabel = 1;
        end
        if strcmp(MODEFL,'classification')
            
            switch dat.TrainParam.PREPROC{1,1}.BINMOD
                case 1
                    inp.multiflag = 0;
                case 0
                    inp.multiflag = 1;
            end
             cv_ExportModels(inp);
        else
            inp.multiflag = 0;
            cv_ExportModels(inp);
        end
    end
    
end
end