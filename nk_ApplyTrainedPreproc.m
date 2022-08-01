function [ inp, contfl, analysis, mapY, GD, MD, Param, P, mapYocv ] = nk_ApplyTrainedPreproc(analysis, inp, paramfl, Param)
% =========================================================================
% [ contfl, analysis, mapY, GD, MD, Param, P, mapYocv ] = ...
%                           nk_ApplyTrainedPreproc2(analysis, inp, paramfl)
% =========================================================================
% Main function to compute /load and return preprocessing parameters and
% preprocessed data for trained analysis chains. The functions is used by
% nk_VisModels and nk_OOCV. It identifies the smallest number of preprocessing 
% parameter combinations needed for a trained model, for which the preprocessing 
% chains have to be computed. This set of parameter combinations is then
% forwarded to nk_PerfPreprocess or nk_Prepwhere the actual computations are
% coordinated.
% 
% Input:
% -------
%
% Output:
% -------
%
% =========================================================================
% (c) Nikolaos Koutsouleris, 12/2021

global VERBOSE PREPROC SAV OCTAVE

contfl = false; 
mapYocv = []; mapY = []; GD = []; MD = []; P = paramfl; 
% Load CVdatamat for current CV2 partition
if ~exist(analysis.RootPath,'dir'), analysis.RootPath = nk_DirSelector('Specify root directory of analysis'); end
if isempty(analysis.GDfilenames{inp.f,inp.d}), contfl = true; return; end
GDpath = fullfile(analysis.RootPath, analysis.GDfilenames{inp.f,inp.d});
if ~exist(GDpath,'file'), fprintf('\n%s not found! Skipping CV2 partition [%g,%g]:\n%s',GDpath, inp.f, inp.d); contfl = true; end

smoothonly = false; if isfield(inp,'smoothonly') && inp.smoothonly, smoothonly = true; end
% Loop through modalities ( if needed )

if isfield(inp,'loadGD') && ~isempty(inp.loadGD) && inp.loadGD
    if VERBOSE, fprintf('\nLoading CVdatamat for CV2 partition [%g,%g]:\n%s', inp.f, inp.d, GDpath); end
    T = load(GDpath);
    if isfield(T,'GD'), GD = T.GD; end
    if isfield(T,'MD'), MD = T.MD; end
end

if smoothonly
    inp = nk_PerfInitSpatial(analysis, inp, paramfl);
    inp.smoothonly = false;
    inp.issmoothed = true;
else
    Yocv = []; 
    issmoothed = false; if isfield(inp,'issmoothed') && inp.issmoothed, issmoothed = true; end
    nM = numel(inp.PREPROC);
    if ~isfield(inp,'saveparam'), inp.saveparam = false; end

    P = cell(1, nM); mapY = cell(1, nM); mapYocv = cell(1,nM); 
    % Check whether optimized preprocessing params exist
    if exist("Param",'var') && ~isempty(Param)
        if iscell(paramfl)
            paramfl{1}.found= true;
        else
            paramfl.found = true;
        end
    else
        Param = cell(1,nM);
        if isfield(inp,'optpreprocmat') && ~isempty(inp.optpreprocmat) && ~inp.saveparam
            if exist(inp.optpreprocmat{inp.f,inp.d},'file')
                fprintf('\nLoading optimized pre-processing parameters for CV2 [%g,%g]:\n%s', ...
                        inp.f, inp.d, inp.optpreprocmat{inp.f,inp.d}); 
                load(inp.optpreprocmat{inp.f,inp.d}); paramfl.found = true;
            else
                if VERBOSE, fprintf('ERROR: Loading of pre-computed parameters not possible because path to file does not anymore exist. Update your paths!'); end
            end
        else
            if VERBOSE, fprintf('\nComputing pre-processing parameters for CV2 [%g,%g].\n', inp.f, inp.d); end
        end
    end

    for n=1:nM
        
        % Check if the modality has to be skipped in case all observations
        % consist of NaN

        % Has the [imaging] data been smoothed already?
        if issmoothed, sstr='s'; else, sstr=''; end

        % Training/CV data extraction from container
        Y = inp.X(n).([ sstr 'Y']);

        % Is the a second Yocv container (e.g. if OOCV data is interpreted, 
        % and multiple OOCV case instances have been generated in nk_MLInterpreter)
        if isfield(inp.X(n),[ sstr 'Yocv2']) && ~isempty(inp.X(n).([ sstr 'Yocv2']))
            Yocv = inp.X(n).([ sstr 'Yocv2']);
        elseif isfield(inp.X(n),[ sstr 'Yocv']) && ~isempty(inp.X(n).([ sstr 'Yocv']))
            Yocv = inp.X(n).([ sstr 'Yocv']); 
        end

        % Do we need to extract a weighting mask, too?
        if isfield(inp.X(n), [ sstr 'Yw']) && ~isempty(inp.X(n).([ sstr 'Yw']))
            inp.iYw = inp.X(n).([ sstr 'Yw']);
        end
    
        % Are there modality-specific runtime parameters that should be extracted?
        if iscell(paramfl)
            tparamfl = paramfl{n};
        else
            tparamfl = paramfl;
        end
        
        % ... and modality-specific preprocessing params?
        if iscell(inp.PREPROC)
            PREPROC = inp.PREPROC{n};
        else
            PREPROC = inp.PREPROC;
        end
        tparamfl.PV = inp.X(n);

        if VERBOSE, fprintf('\nGenerate pre-processing parameter array for CV2 partition [%g,%g].\n',inp.f,inp.d); end
        tparamfl = nk_PrepPreprocParams(PREPROC, tparamfl, analysis, n, inp.ll, inp.curlabel);

        % Param is a structure that contains all relevant info to generate the features 
        % needed by the optimized classifier / predictor system
        if tparamfl.found 
            tparamfl.Param = Param{n};
        elseif isfield(tparamfl,'Param')
            tparamfl = rmfield(tparamfl,'Param');
        end
        
        % now we are ready to run the preprocessing pipeline
        if inp.stacking
            [mapY{n}, Param{n}, P{n}, mapYocv{n}] = nk_PerfPreprocessMeta(inp, inp.label, tparamfl);
        else
            [mapY{n}, Param{n}, P{n}, mapYocv{n}] = nk_PerfPreprocess(Y, inp, inp.label, tparamfl, Yocv);
        end

    end
        
    % Save parameters to disk
    if isfield(inp,'saveparam') && inp.saveparam
        CV2p = inp.f; CV2f = inp.d; 
        if isfield(inp,'CV1')
            CV1p = inp.CV1p;  CV1f = inp.CV1f;
        else
            CV1p = []; CV1f = [];
        end
        OptPreprocParamFilename = ...
            nk_GenerateNMFilePath( inp.saveoptdir, SAV.matname, 'OptPreprocParam', [], inp.varstr, inp.id, CV2p , CV2f, CV1p, CV1f );
        fprintf('\nSaving %s to disk...', OptPreprocParamFilename)
        if ~exist(inp.saveoptdir,'dir'), mkdir(inp.saveoptdir); end
        ofold = CV2f; operm = CV2p;
        if OCTAVE
            save(OptPreprocParamFilename,'Param','ofold','operm');
        else
            save(OptPreprocParamFilename,'-v7.3', 'Param','ofold','operm');     
        end
    end

    % Transfer mapped data to appropriate container
    if nM > 1
       mapY = nk_mapY2Struct(mapY, true);
       if (iscell(mapYocv) && ~sum(cellfun(@isempty,mapYocv))) || ( ~iscell(mapYocv) && ~isempty(mapYocv))
           mapYocv = nk_mapY2Struct(mapYocv, true); 
       end
    else
        mapY = mapY{n};
        if ~isempty(mapYocv), mapYocv = mapYocv{n}; end
    end
end