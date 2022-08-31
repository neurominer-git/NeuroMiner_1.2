function [act, SEQOPT ] = nk_SEQOPT_config(SEQOPT, defaultsfl)

if ~exist('defaultsfl','var') || isempty(defaultsfl), defaultsfl = false; end

if ~defaultsfl 
    
    if ~isfield(SEQOPT,'C'),            SEQOPT.C = NaN; end
    if ~isfield(SEQOPT,'PerfMode'),     SEQOPT.PerfMode = 1; end
    if ~isfield(SEQOPT,'AnchorType'),   SEQOPT.AnchorType = 1; end
    if ~isfield(SEQOPT,'ReplaceMode'),  SEQOPT.ReplaceMode = 1; end
    if ~isfield(SEQOPT,'Mode'),         SEQOPT.Mode = 1; end
    
    MODEARR = {'The sample to be deferred','The entire population'}; MODESTR = MODEARR{SEQOPT.Mode};
    MnuStr = sprintf('Define the population the deferral optimization will be based upon [ %s ]', MODESTR);         MnuAct = 1;
    
    REPLACEARR = {'Sequential (=replacing former with latter predictions in the sequence)','Integrative (=computing the mean across deferral nodes)'};  REPLACESTR = REPLACEARR{SEQOPT.ReplaceMode};
    MnuStr = sprintf('%s|Define how to combine sequential predictions [ %s ]',MnuStr, REPLACESTR); MnuAct = [MnuAct 2];
    
    ANCHORARR = {'Decision boundary','Median of decision score distribution'}; ANCHORSTR = ANCHORARR{SEQOPT.AnchorType};
    MnuStr = sprintf('%s|Define how of the deferral window is anchored [ %s ]',MnuStr, ANCHORSTR);                  MnuAct = [MnuAct 3];
    
    
    [~ , ~, ~, critstr ] = nk_GetScaleYAxisLabel;
    OPTIMARR = {sprintf('Current optimization criterion (%s)',critstr),'Mean decision distance between classes'}; OPTIMSTR = OPTIMARR{SEQOPT.PerfMode};
    MnuStr = sprintf('%s|Define the optimisation criterion for deferral learning [ %s ]', MnuStr, OPTIMSTR);        MnuAct = [MnuAct 4];
    
    if isnan(SEQOPT.C)
        SEQSTR = 'Deferral sequence undefined';
    else
        SEQSTR = sprintf('%g deferral sequence(s) defined, max sequence length: %g predictors',size(SEQOPT.C,1), size(SEQOPT.C,2));
    end
    MnuStr = sprintf('%s|Define the deferral sequences to be tested [ %s ]', MnuStr, SEQSTR);                       MnuAct = [MnuAct 5];
    
    nk_PrintLogo
    act = nk_input('Select parameters for deferral sequence optimization',0,'mq', MnuStr, MnuAct, 1);
    switch act
        case 1
            SEQOPT.Mode         = nk_input('Define training population for deferral optimization',0,'m','The sample to be deferred|The entire population',[1,2],SEQOPT.Mode);
        case 2
            SEQOPT.ReplaceMode  = nk_input('Sequential decision computation',0,'m','Sequential|Integrative',[1,2],SEQOPT.ReplaceMode);
        case 3
            SEQOPT.AnchorType   = nk_input('Anchor deferral window to model''s decision boundary or decision score median',0,'m','Decision boundary|Median',[1,2], SEQOPT.AnchorType);
        case 4
            CurrOptStr = sprintf('Current optimization criterion (%s)', critstr);
            SEQOPT.PerfMode     = nk_input('Deferral optimization mode',0,'m',[ CurrOptStr '|Mean decision distance'],[1,2], SEQOPT.PerfMode);
        case 5
            SEQOPT.C            = nk_input('Define matrix of deferral sequences to be evaluated',0,'e');
    end
else
    act =0; SEQOPT = struct('C', NaN, 'PerfMode', 1, 'AnchorType', 1, 'ReplaceMode', 1, 'Mode', 1); 
end


