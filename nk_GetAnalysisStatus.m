function Status = nk_GetAnalysisStatus(NM)

Status.completed_analyses  = [];  
Status.isequal_cv          = []; 
Status.nmodal_analyses     = [];
Status.analexistflag       = false; 
Status.analreadyflag       = false; 
Status.analcompleteflag    = false;
Status.oocvreadyflag       = false;
Status.oocvappflag         = false;
Status.stacking_analyses   = [];
Status.n_inputanalyses     = [];
Status.sequence_analyses   = [];

if isfield(NM,'analysis')
    n_anal = numel(NM.analysis);
    Status.completed_analyses = false(1,n_anal);
    Status.isequal_cv = false(1,n_anal);
    Status.nmodal_analyses = zeros(1,n_anal);
    Status.sequence_analyses = false(1,n_anal);
    Status.stacking_analyses = false(1,n_anal);
    Status.n_inputanalyses = zeros(1,n_anal);
    for i = 1:n_anal
        if NM.analysis{i}.status
            Status.completed_analyses(i)= true; 
        end
        if isequaln(NM.cv, NM.analysis{i}.params.cv)
            Status.isequal_cv(i) = true; 
        end
        if Status.completed_analyses(i)
            Status.nmodal_analyses(i) = numel(NM.analysis{i}.GDdims);
        end
        if NM.analysis{i}.params.TrainParam.STACKING.flag==1
            Status.stacking_analyses(i) = true;
            Status.n_inputanalyses(i) = numel(NM.analysis{i}.params.TrainParam.STACKING.sel_anal); 
        end
        if strcmp(NM.analysis{i}.params.TrainParam.SVM.prog,'SEQOPT')
            Status.sequence_analyses(i) = true;
        end
    end
end

if ~isempty(Status.completed_analyses)    
    Status.analexistflag = true; 
    if any(Status.completed_analyses), Status.analreadyflag = true; end    
    if sum(Status.completed_analyses) == numel(Status.completed_analyses), Status.analcompleteflag = true; end
    if Status.analcompleteflag && isfield(NM,'OOCV'), Status.oocvreadyflag = true; end
    if Status.oocvreadyflag && isfield(NM.defs,'analyses_locked') && NM.defs.analyses_locked, Status.oocvappflag = true; end
end

