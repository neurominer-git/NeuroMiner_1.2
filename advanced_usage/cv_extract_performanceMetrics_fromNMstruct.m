function perftab = cv_extract_performanceMetrics_fromNMstruct(NMstruct, varargin)
if nargin>1
   
    analyses_status = nk_GetAnalysisStatus(NMstruct);
    analyses_ind = find(analyses_status.completed_analyses);
    if sum(analyses_ind) == 0
        warning('no analyses completed'); 
        exit; 
    end
    oocv_flag = 0; 
    save_flag = 0; 
    metric_selection = fieldnames(NMstruct.analysis{1, analyses_ind(1)}.GDdims{1, 1}.BinClass{1, 1}.contigency);
    i = 1;
    while i < length(varargin)
        if ischar(varargin{i})
            switch varargin{i}
                case 'analyses_ind'
                    suggested_analyses_ind = varargin{i+1};
                    if all(ismember(suggested_analyses_ind, analyses_ind))
                        analyses_ind = suggested_analyses_ind; 
                    end
                    i = i+2; 
                case 'oocv_container'
                    oocv_flag = 1;      
                    oocv_container = varargin{i+1};
                    i = i+2;
                case 'metric_selection'
                    suggested_metric_selection = varargin{i+1}; 
                    if all(ismember(suggested_metric_selection, metric_selection))
                        metric_selection = suggested_metric_selection; 
                    end
                    i = i+2; 
                case 'filename'
                    save_flag = 1; 
                    filename = varargin{i+1};
                    i = i+2; 
                otherwise 
                    warning('unknown option'); 
                    i = i+2; 
            end
        else 
            warning('wrong input argument'); 
        end
    end


    perftab = array2table(zeros(length(analyses_ind),length(metric_selection)), 'VariableNames', metric_selection);
    analysis_IDs = [];
    if oocv_flag
        for a = 1:length(analyses_ind)
            analysis = analyses_ind(a);
            analysis_IDs = [analysis_IDs, {NMstruct.analysis{1,analysis}.id}];
            for m = 1:length(metric_selection)
                metric = metric_selection{m};
                eval(sprintf("perftab{a,metric} = NMstruct.analysis{1, analysis}.OOCV{1, oocv_container}.BinResults{1, 1}.contingency{1, 1}.%s", metric));
            end
        end
    else
        for a = 1:length(analyses_ind)
            analysis = analyses_ind(a);
            analysis_IDs = [analysis_IDs, {NMstruct.analysis{1,analysis}.id}];
            for m = 1:length(metric_selection)
                metric = metric_selection{m};
                eval(sprintf("perftab{a,metric} = NMstruct.analysis{1, analysis}.GDdims{1,1}.BinClass{1, 1}.contigency.%s", metric));
            end
        end
    end
    perftab.Analysis_ID = analysis_IDs'; 

    
    if save_flag
        writetable(perftab, filename, "Delimiter",",");
    end

end