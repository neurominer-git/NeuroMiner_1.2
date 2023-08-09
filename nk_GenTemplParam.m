function TEMPL = nk_GenTemplParam(PREPROC, CV, MODEFL, RAND, Y, inp, kbin, paramfl, act)

if ~exist("act","var") || isempty(act), act = 'params'; end

% For factorization methods: TEMPLATE MAPPING 
if PREPROC.BINMOD
    ukbin = kbin;   SrcParam.binmult = 1;
else
    ukbin = 1;      SrcParam.binmult = 0;
end
SrcParam.CV1perm            = 1; 
SrcParam.CV1fold            = 1;
SrcParam.covars             = inp.covars;
SrcParam.covars_oocv        = inp.covars_oocv;

fprintf('\nGenerating template preprocessing parameters using entire dataset.')

for curclass = 1 : ukbin

    SrcParam.u = curclass;
    paramfl.P{curclass} = nk_ReturnParamChain(PREPROC, 1);
    if ~ukbin>1 || strcmp(MODEFL,'regression')
        SrcParam.TrX = true(length(inp.labels),1);
    else
        SrcParam.TrX = find(inp.label == CV.class{1,1}{curclass}.groups(1) | inp.label == CV.class{1,1}{curclass}.groups(2));   
    end
    if iscell(Y)
        iY = Y{curclass}{1}(SrcParam.TrX,:);
    else
        iY = Y(SrcParam.TrX,:);
    end
    [InputParam.Tr , TrLX, SrcParam.iTrX] = nk_ManageNanCases(iY, SrcParam.TrX); 
    switch MODEFL
        case 'classification' 
             if RAND.Decompose ~=9
                SrcParam.BinaryTrainLabel   = SrcParam.TrX;
                SrcParam.BinaryCVLabel      = SrcParam.TrX;
            end
            SrcParam.MultiTrainLabel    = inp.label;
            SrcParam.MultiCVLabel       = inp.label;
        case 'regression'
            SrcParam.TrainLabel         = inp.label;
            SrcParam.CVLabel            = inp.label;
    end
    if isfield(inp,'Yw'), InputParam.Yw = inp.Yw; end 
    [TEMPL.Tr{curclass}, oTrainedParam] = nk_GenPreprocSequence(InputParam, PREPROC, SrcParam);

    if isfield(paramfl,'PREPROC') && isfield(paramfl,'PXfull') && ~isempty(paramfl.PXopt{curclass})
        % Here an optimized parameter space exists that has
        % been used to limit the computation to the unique
        % parameter cominations. We create and store the pointers
        % and used them later on to retrieve the right preproc
        % version of the data and the respective preproc params
        [ TEMPL.Param(1,1,curclass).data_ind, ...
          TEMPL.Param(1,1,curclass).train_ind, ...
          TEMPL.Param(1,1,curclass).nP, ...
          TEMPL.Param(1,1,curclass).nA] = nk_ParamReplicator(paramfl.P{curclass}, paramfl.PXopt{curclass}, paramfl.PREPROC, numel(oTrainedParam));
    end
    TEMPL.Param(1,1,curclass).TrainedParam = oTrainedParam;

end
if strcmp(act,'params'), TEMPL.Tr = []; end