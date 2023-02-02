function TEMPL = nk_GenTemplParam(PREPROC, CV, MODEFL, RAND, Y, inp, kbin)

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
   
    switch MODEFL
        case 'classification' 
            SrcParam.TrX = find(inp.label == CV.class{1,1}{curclass}.groups(1) | inp.label == CV.class{1,1}{curclass}.groups(2));
            if RAND.Decompose ~=9
                SrcParam.BinaryTrainLabel   = SrcParam.TrX;
                SrcParam.BinaryCVLabel      = SrcParam.TrX;
            end
            SrcParam.MultiTrainLabel    = inp.label;
            SrcParam.MultiCVLabel       = inp.label;
        case 'regression'
            SrcParam.TrX = true(length(inp.labels),1);
            SrcParam.TrainLabel         = inp.label;
            SrcParam.CVLabel            = inp.label;
    end
     if iscell(Y)
        iY = Y{curclass}{1}(SrcParam.TrX,:);
    else
        iY = Y;
    end
    [InputParam.Tr , SrcParam.TrX, SrcParam.iTr] = nk_ManageNanCases(iY, SrcParam.TrX); 
    if isfield(inp,'Yw'), InputParam.Yw = inp.Yw; end 
    [TEMPL.Tr{curclass}, TEMPL.Param{curclass}] = nk_GenPreprocSequence(InputParam, PREPROC, SrcParam);

end