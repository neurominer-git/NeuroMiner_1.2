function [ GD, MD ] = nk_MLOptimizer_ParamCycler(GD, MD, DISP, Ps, Params_desc, mapY, algostr, f, d, npreml, nclass, ngroups, batchflag, PsSel, combcell)
% =========================================================================
% FORMAT [ GD, MD ] = nk_MLOptimizer_ParamCycler(GD, MD, DISP, Ps, ...
%                           Params_desc, mapY, algostr, f, d, npreml, ...
%                           nclass, batchflag, PsSel, combcell)
% =========================================================================
% This child function of nk_MLOptimizer performs an extensive (brute-force)
% search of the parameter space as defined by the user.
% 
% Inputs:
% -------
% GD            : Results container
% MD            : Model container
% DISP          : Display structure with data for NM Optimization Viewer
% Ps            : Parameter combinations
% Params_desc   : Descriptions of parameters
% mapY          : The data containing CV1 training and CV1 test and CV2
%                 validation data
% algostr       :
% [ f, d ]      : Position in the CV2 grid
% npreml        : Number of free preprocessing parameters in the space
% nclass        : Number of binary classifier | predictors
% batchflag     : batchmode (for HPC) ?
% PsSel         : Previously selected parameters nodes and data at that
%                 nodes
% combcell      : Flag indicating that Ps is a cell array rather than a
%                 numeric array of free parameters
%
% Outputs:
% --------
% GD (see above)
% MD (see above)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Nikolaos Koutsouleris, 11/2022

global CV MULTILABEL CVPOS

nPs = size(Ps{1},1); 
CVPOS.CV2p = f;
CVPOS.CV2f = d;

% Any free preprocessing parameters ?
if npreml>-1
    if combcell
        pp = unique(cell2mat(Ps{1}(:,end-npreml:end)),'rows','stable');
    else
        pp = unique(Ps{1}(:,end-npreml:end),'rows','stable');
    end
end

% Do we have a multi-label situation?
nl = nk_GetLabelDim(MULTILABEL);

tElapsedSum = 0; ii = [];
tic
labelstr = ''; 
if nPs>1, fprintf('\n === Performing hyperparameter optimization === \n'); else, fprintf('\n'); end
     
% Loop through all available labels
for curlabel=1:nl
    if MULTILABEL.flag
        if nl>1
            labelstr = sprintf('Label #%g: %s | ', curlabel, MULTILABEL.desc{curlabel});
        else
            labelstr = sprintf('Label %s | ', MULTILABEL.desc{curlabel});
        end
    end
    MULTILABEL.curdim = curlabel;
    
    % Check whether selected parameter combinations should be tested
    if ~exist('PsSel','var') || isempty(PsSel)
        PiSel = true(nPs,nclass);
    else
        PiSel = false(nPs,nclass);
        for curclass=1:nclass
            PiSel(:,curclass) = PsSel{curclass}{curlabel}.SelNodes;
        end 
    end
    pltmax = sum(any(PiSel,2));
    
    % Do we have different CV structures assigned to different labels
    TCV = CV;
    if numel(TCV)>1, CV = TCV(curlabel); end
    pltcnt =0 ; 
    
    % Loop through all parameter combinations
    for i = 1:nPs

        if ~sum(any(PiSel(i,:)))
            sskip = sprintf('\n%sSkipping hyperparameter combination %g!',labelstr,i);
            fprintf('%s',sskip); %fprintf(repmat('\b',1,numel(sskip))); 
            continue; 
        end
        ii = [ii i];
        pltcnt = pltcnt+1; pltperc = pltcnt*100/pltmax;
        tElapsed = toc; tElapsedSum = tElapsedSum+tElapsed; 
        elaps = sprintf('\t%1.2f sec.',tElapsed);
        if nPs > 1 
            DISP.s = sprintf('%s | %s%s\nCV2 [ %g, %g ] => %4g/%4g parameter combinations => %1.1f%% ', ...
                elaps, labelstr, algostr, f, d , pltcnt, pltmax, pltperc);
        else
            DISP.s = sprintf('%s | %s%s\nCV2 [ %g, %g ] => No-parameter optimization', ...
                elaps, labelstr, algostr, f, d  );
        end
        fprintf('%s',DISP.s); 
        DISP.pltperc = pltperc; 

        %% Loop through binary learners and prepare learning params
        tic;
        cPs = cell(nclass,1);
        for curclass = 1:nclass
            DISP.P{curclass} = Ps{curclass}(i,:);
            cPs{curclass} = nk_PrepMLParams(Ps{curclass}, Params_desc{curclass}, i);
        end

        %% Check whether new mapY container has to be retrieved
        % Now retrieve preprocessing parameter combinations and check whether
        % there has been a parameter change from the previous to current
        % parameter combination
        dimchng = false; i_dl = 1; m_dl = 1;
        if npreml >-1
            if numel(ii) > 1
                if combcell
                    ipp     = cell2mat(Ps{1}(ii(end),end-npreml:end));
                    impp    = cell2mat(Ps{1}(ii(end-1),end-npreml:end));
                else
                    ipp     = Ps{1}(ii(end),end-npreml:end);
                    impp    = Ps{1}(ii(end-1),end-npreml:end);
                end
                [~,i_dl] = ismember(ipp,pp,'rows');
                [~,m_dl] = ismember(impp,pp,'rows');
            else
                if combcell
                    ipp     = cell2mat(Ps{1}(ii,end-npreml:end));
                else
                    ipp     = Ps{1}(ii,end-npreml:end);
                end
                [~,i_dl] = ismember(ipp,pp,'rows');
                m_dl = i_dl;
            end
        end
        if i_dl ~= m_dl || (~exist('mapYi','var') || isempty(mapYi))
            dimchng = true;
        end

        %% Model training phase
        if dimchng % now retrieve new mapYi from container
            mapYi = nk_MLOptimizer_ExtractDimMat(mapY, i_dl, cPs); 
            % ... and extract features according to filter mechanism (if
            % needed)
            FilterSubSets = nk_CreateSubSets(mapYi); 
        end   
        % Compute current model(s) for variable parameter combination P(i) = [ P1 P2 P3
        % ... Pn] using the CV1 data partitions. Apply single or ensemble model
        % to CV2 test sample in order to estimate the generalization 
        % capacity of the classification / prediction rule    
        [CV1perf, CV2perf, models] = nk_CVPermFold(mapYi, nclass, ngroups, cPs, FilterSubSets, batchflag);      

        % Transfer results from CV1perf and CV2perf to GD
        % structure using nk_GridSearchHelper2 function
        [GD, MD, DISP] = nk_GridSearchHelper(GD, MD, DISP, i, nclass, ngroups, CV1perf, CV2perf, models);
       
        if isfield(CV1perf,'detrend'), GD.Detrend{i} = CV1perf.detrend; end

        % Create variate mask according to selected features
        if isfield(mapYi,'VI')

            [iy,jy] = size(CV(1).cvin{f,d}.TrainInd);
            
            GD.VI{i,curlabel} = cell(iy,jy,nclass);
            for k=1:iy
                for l=1:jy
                    if iscell(mapYi.VI{k,l})
                         for curclass = 1:nclass
                            VI = repmat(mapYi.VI{k,l}{curclass},1,size(GD.FEAT{i,curlabel}{k,l,curclass},2));
                            GD.VI{i,curlabel}{k,l,curclass} = VI;
                         end
                    else
                        for curclass = 1:nclass
                            VI = repmat(mapYi.VI{k,l},1,size(GD.FEAT{i,curlabel}{k,l,curclass},2));
                            GD.VI{i,curlabel}{k,l,curclass} = VI;
                        end
                    end                            
                end
            end
            
        end
        fprintf(repmat('\b',1,numel(DISP.s))); 
    end
end
MULTILABEL.curdim = 1;
CV = TCV(1);

fprintf('\n');fprintf('CV2 [%g, %g]: OPTIMIZATION COMPLETED IN %1.2f SEC ', ...
    f, d, tElapsedSum)
