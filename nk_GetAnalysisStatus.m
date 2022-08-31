function Status = nk_GetAnalysisStatus(NM)

Status.completed_analyses       = [];  
Status.isequal_cv               = []; 
Status.nmodal_analyses          = [];
Status.analexistflag            = false; 
Status.analreadyflag            = false; 
Status.analcompleteflag         = false;
Status.oocvreadyflag            = false;
Status.oocvappflag              = false;
Status.stacking_analyses        = [];
Status.n_inputanalyses          = [];
Status.sequence_analyses        = [];
Status.analyses_nondeterministic = [];
Status.analyses_visualized      = [];
Status.analyses_imaging         = [];
Status.analyses_interpreted     = [];

if isfield(NM,'analysis')
    
    n_anal = numel(NM.analysis);
    Status.completed_analyses       = false(1,n_anal);
    Status.isequal_cv               = false(1,n_anal);
    Status.nmodal_analyses          = zeros(1,n_anal);
    Status.sequence_analyses        = false(1,n_anal);
    Status.stacking_analyses        = false(1,n_anal);
    Status.n_inputanalyses          = zeros(1,n_anal);
    Status.analyses_visualized      = false(1,n_anal);
    Status.analyses_imaging         = cell(1,n_anal);
    Status.analyses_interpreted     = false(1,n_anal);

    for i = 1:n_anal
        % Completed analyses
        if NM.analysis{i}.status
            Status.completed_analyses(i)= true; 
        end
        % Analyses with equal CV structures (according to NM workspace CV)
        if isequaln(NM.cv, NM.analysis{i}.params.cv)
            Status.isequal_cv(i) = true; 
        end
        % Number of modalities
        if Status.completed_analyses(i)
            Status.nmodal_analyses(i) = numel(NM.analysis{i}.GDdims);
        end
        % Stacking analyses
        if NM.analysis{i}.params.TrainParam.STACKING.flag==1
            Status.stacking_analyses(i) = true;
            Status.n_inputanalyses(i) = numel(NM.analysis{i}.params.TrainParam.STACKING.sel_anal); 
        end
        % Sequence analyses
        if strcmp(NM.analysis{i}.params.TrainParam.SVM.prog,'SEQOPT')
            Status.sequence_analyses(i) = true;
        end
        % Fusion analyses
        FUSION = NM.analysis{i}.params.TrainParam.FUSION;
        Modality = FUSION.M;
        switch FUSION.flag
            case {0, 2, 3}
                nM = numel(Modality);
            case 1
                nM = 1;
        end
        % Non-deterministic analyses
        if isfield(NM.analysis{i}.params.TrainParam.SVM,'ADASYN') && NM.analysis{i}.params.TrainParam.SVM.ADASYN.flag==1
            for j=1:nM
                PREPROC = NM.analysis{i}.params.TrainParam.PREPROC{Modality(j)};
                Status.analyses_nondeterministic(i).modality(j) = false;
                if isfield(PREPROC,'ACTPARAM')
                    for q=numel(PREPROC.ACTPARAM):-1:1
                        if strcmp(PREPROC.ACTPARAM{q}.cmd,'reducedim') 
                            if strcmp(PREPROC.ACTPARAM{q}.DR.RedMode,'PCA') && PREPROC.ACTPARAM{q}.DR.PercMode == 3
                                Status.analyses_nondeterministic(i).modality(j) = true;
                            end
                        end
                    end
                end
            end
        else
           for j=1:nM
               Status.analyses_nondeterministic(i).modality = false(1,nM);
           end
        end
        if isfield(NM.analysis{i},'visdata')
            Status.analyses_visualized(i) = true;
        end
        for j=1:nM
            if strcmp(NM.analysis{i}.params.datadescriptor{Modality(j)}.source,'image')
                Status.analyses_imaging{i} = [ Status.analyses_imaging{i}  Modality(j)];
            else
                Status.analyses_imaging{i} = [];
            end
        end
        if isfield(NM.analysis{i},'MLI')
            Status.analyses_interpreted(i) = true;
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

