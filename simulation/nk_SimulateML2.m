function [ Res, IN ] = nk_SimulateML2(IN, axes)
% =========================================================================
% FORMAT function [ R, P ] = nk_SimulateML(IN)
% =========================================================================
% This function is part of the NeuroMiner machine learning library.
% It allows the user to determine the expected prognostic performance of a
% binary classification algorithm under different learning scenarios, as
% determined by the "IN" parameter structure:
%
% 1) IN.Nfeats : the dimensionality of the variable/feature space
%
% 2) IN.Nmarkers : the percentage (0-1) of predictive markers within that
%    feature space
%
% 3) IN.Ncases : the number of observations/cases/samples to learn from
%
% 4) IN.eventprob : the percentage (0-1) of cases with the (desired) event
%    in relation to the total number of cases
%
% 5) IN.AUCmax and IN.AUCmin : the maximum / minimum separability in the
%    feature matrix. The helper function is nk_ExpVec.m is used to create
%    features following an exponntial decline in separability. Feature
%    creation is based on randomly picking Gaussians within a certain
%    distance of each other and a certain span.
%
% 6) IN.Nbatches : the absolute number of batched to corrupt the data
%
% 7) IN.BatchPerc : the maximum percentage (0-1) of the feature value range
%    allowed for site effects
%
% 8) IN.NcasesMiss : The percentage (0-1) of cases with missing data
%
% 9) IN.NfeatsMiss : The percentage (0-1) of missing values per case
%
% 10) IN.Data : The original data which structure should be simulated
%
% 11) IN.DataLabel : Output variable of original data
%
% 12) IN.NRanalysis : Which analysis to simulate, index (from NM structure)
%
% Parameters can come in ranges and thus are hyperparameters that the
% function will loop through in a brute-force approach. Based on these
% hyperparameters the function will generate a synthetic data matrix and
% forward it to a simple machine learning pipeline. The user currently can
% choose among four different classifiers:
% IN.algorithm :
%       a) "LINKERNSVM" : LIBSVM with linear kernel (C=1)
%       b) "LINSVM": LINBLINEAR, L2-regularized, L1-Loss SVC
%       c) "L2LR" : L2-regularized logistic regression
%       d) "L1LR" : L1-regularized logistic regression
%       [not implemented yet:  "RF" : Random Forests]
%
% Furthermore, the setup of the cross-validation cycle can be tweaked, e.g.:
% IN.RAND.OuterPerm = 1; ==> Outer cross-validation permutations
% IN.RAND.InnerPerm = 1; ==> Inner cross-validation permutations
% IN.RAND.OuterFold = 5; ==> Outer cross-validation folds
% IN.RAND.InnerFold = 5; ==> Inner cross-validation folds
% IN.RAND.Decompose = 1;
%
% Please note that currently the function does not run any nested
% cross-validation as we will use standard machine learning parameters.
% This could be changed in the future but it will increase computational
% costs of the simulation enormously.
%
% This function will only work with NeuroMiner installed and initialized in
% an active MATLAB session Use preferably MATLAB 2018b or higher.
% _________________________________________________________________________
% (c) Nikolaos Koutsouleris, 10/02/2022
global fromData xNM

if isfield(IN,'Data') && ~isempty(IN.Data)
    fromData = 1;
    xNM = evalin('base','NM');
    analysis = xNM.analysis{1,IN.NRanalysis};
    %     xNM.analysis = cell(1,size(IN.Ncases,2));
    %     xNM.Y = cell(1,size(IN.Ncases,2))
    %     for i = 1:size(IN.Ncases,2)
    %         xNM.Y{1,i} =
    %     end
else
    fromData = 0;
end

% Create hyperparameter array
% Default algorithm
if ~isfield(IN,'algorithm') || isempty(IN.algorithm)
    IN.algorithm = 'LINKERNSVM';
end

if ~isfield(IN,'Nbatches') || isempty(IN.Nbatches)
    IN.Nbatches = 0;
end

if ~isfield(IN,'BatchPerc') || isempty(IN.BatchPerc)
    IN.BatchPerc = 0.1;
end

if ~isfield(IN,'NcasesMiss') || isempty(IN.NcasesMiss)
    IN.NcasesMiss = 0;
end

if ~isfield(IN, 'NfeatsMiss') || isempty(IN.NfeatsMiss)
    IN.NfeatsMiss = 0;
end

if ~isfield(IN, 'AUCmax') || isempty(IN.AUCmax)
    IN.AUCmax = 0.7;
end

if ~isfield(IN, 'AUCmin') || isempty(IN.AUCmin)
    IN.AUCmin = 0.5;
end

if ~isfield(IN, 'verbose') || isempty(IN.verbose)
    IN.verbose = true;
end

if fromData
    origRAND = analysis.params.TrainParam.RAND;
    P = IN.Ncases;
    nP = size(P,2);
else
    if ~isfield(IN, 'RAND') || isempty(IN.RAND)
        % Setup permutation parameters
        IN.RAND.OuterPerm = 1;
        IN.RAND.InnerPerm = 1;
        IN.RAND.OuterFold = 5;
        IN.RAND.InnerFold = 5;
        IN.RAND.Decompose = 1;

    end
    P = allcomb2(IN.Nfeats, IN.Nmarkers, IN.Ncases, IN.eventprob, IN.AUCmax, IN.AUCmin, IN.Nbatches, IN.BatchPerc, IN.NcasesMiss, IN.NfeatsMiss);
    nP = size(P,1);
end


R = zeros(nP,1);
R95CI = zeros(2,nP);


if nargin == 2
    ax = axes;
else
    figure;
    ax = gca;
end
Crit = 'Balanced Accuracy';
Xl = cell(nP,1);


for i=1:nP
    if fromData
        sC = sprintf('C: %g/%g', i, size(P,1));
    else
        if numel(IN.Nfeats)>1,       sF = sprintf('F%g', P(i,1)); else, sF= ''; end
        if numel(IN.Nmarkers)>1,     sM = sprintf(' M:%g', P(i,2)); else, sM= ''; end
        if numel(IN.Ncases)>1,       sC = sprintf(' C:%g', P(i,3)); else, sC= ''; end
        if numel(IN.eventprob)>1,    sE = sprintf(' E:%g', P(i,4)); else, sE= ''; end
        if numel(IN.AUCmax)>1,       sAU = sprintf(' AU:%g', P(i,5)); else, sAU= ''; end
        if numel(IN.AUCmin)>1,       sAL = sprintf(' AL:%g', P(i,6)); else, sAL= ''; end
        if numel(IN.Nbatches)>1,     sB = sprintf(' B:%g', P(i,7)); else, sB= ''; end
        if numel(IN.BatchPerc)>1,    sBP = sprintf(' BP:%g', P(i,8)); else, sBP= ''; end
        if numel(IN.NcasesMiss)>1,   sCm = sprintf(' Cm:%g', P(i,9)); else, sCm= ''; end
        if numel(IN.NfeatsMiss)>1,   sFm = sprintf(' Fm:%g', P(i,10)); else, sFm= ''; end
        Xl{i} = sprintf('%s%s%s%s%s%s%s%s%s%s',sF, sM, sC, sE, sAU, sAL, sB, sBP, sCm, sFm);
    end
end
for i=1:nP % Loop through hyperparameter combinations
    fprintf('\nWorking on parameter combination %g/%g: %s', i, nP, Xl{i});
    % Do the magic
    if fromData
        [R(i), R95CI(:,i)] = compute_perf(IN.Nfeats, ...
            IN.Nmarkers, ...
            IN.Ncases(:,i), ...
            IN.eventprob, ...
            IN.AUCmax, ...
            IN.AUCmin, ...
            IN.Nbatches, ...
            IN.BatchPerc, ...
            IN.NcasesMiss, ...
            IN.NfeatsMiss, ...
            IN.algorithm, ...
            origRAND, ...
            IN.verbose, IN.Modalities, IN.DataLabel, IN.NRanalysis, IN.add2orig, IN.NReps, IN.SitesIdx, IN.condIdx); % varargin
    else
        [R(i), R95CI(:,i)] = compute_perf(P(i,1), P(i,2), P(i,3), P(i,4), P(i,5), P(i,6), P(i,7), P(i,8), P(i,9), P(1,10), IN.algorithm, IN.RAND, IN.verbose);
    end
    % Update the simulation plot
    plot(ax,1:i,R(1:i),'b-');
    ax.XTick = 1:i;
    ax.XTickLabel = Xl(1:i);
    ax.XTickLabelRotation = 90;
    ax.XAxis.Label.String = 'Parameter combinations'; ax.XAxis.Label.FontWeight='bold';
    ax.YAxis.Label.String = Crit; ax.YAxis.Label.FontWeight='bold';
    drawnow
end
Res.R = R;
Res.R95CI = R95CI';
Res.Params = P;

function [Rmean, R95CI] = compute_perf(nf, mr, nc, er, auc_max, auc_min, nb, bp, ncm, nfm, algorithm, RAND, varargin)
global SVM fromData xNM xCV



% Create marker matrix based on user's input
% fill the rest of the matrix with uninformative features
% in the range of -1 to 1

% check if original data was provided
if fromData
    xNM = evalin('base','NM');
    verbose = varargin{1,1};
    %Ys = varargin{1,2};
    mods = varargin{1,2};
    labels = varargin{1,3};
    curanal = varargin{1,4};
    add2orig = varargin{1,5};
    reps = varargin{1,6};
    sitesIdx = varargin{1,7};
    condIdx = varargin{1,8};
    if condIdx == 0
        condName = 'label';    
    elseif condIdx > 0
        condName = app.NM.covnames(condIdx);
    end

    origRAND = RAND;
    

    % what happens in the next lines, really necessary?
%     if isa(labels_a, 'cell') % only binary problems
%         num_labels = double(strcmp(labels_a,labels_a(1,1)));
%         num_labels(num_labels==0) = -1;
%     else % data and labels from NM structure; binary classification
%         num_labels = labels_a;
%         num_labels= num_labels-1; % simulate function could potentially better be adjusted to fit NM strucutre of labels
%     end


    analrootdir = xNM.analysis{1,curanal}.rootdir;
    origDataFile = sprintf('%s/origData.csv',analrootdir);

    
    % check if several modalities included
    %   if length(Ys) >1
    % mods = 1:length(Ys);
    modNfeats = zeros(1,length(mods));
    Y = [];
    YColNames = [];
    for i = 1:length(mods)
        modNfeats(i) = size(xNM.Y{1,mods(i)},2);
        Y = [Y, xNM.Y{1,i}];   
    end
    modColIdx = repelem(mods, modNfeats);
    modNameVec = "Ymod" + modColIdx + "_Feat" + (1:sum(modNfeats));
    YColNames = [YColNames, modNameVec];
    %   end

    %     if isa(varargin{1,2},'double')
    %         Ytab = array2table(varargin{1,2});
    %     end

    %     Y = varargin{1,2};


    % check whether covariates are included in analysis
    if isfield(xNM, 'covars')
        covColIdx = zeros(1,size(Y,2));
        sitesIdxY = sitesIdx + size(Y,2);
        Y = [Y, xNM.covars];
        covColIdx = [covColIdx, ones(1,size(xNM.covars,2))];
        YColNames = [YColNames, xNM.covnames'];
    end
  
    % check whether leave-one-group out cv framework is selected
    cv1lcoIdx = 0;
    if isfield(origRAND,'CV1LCO')
        cv1lcoColIdx = zeros(1,size(Y,2));
        Y = [Y, origRAND.CV1LCO.ind]
        cv1lcoColIdx =[cv1lcoColIdx, 1];
        cv1lcoIdx = length(cv1lcoColIdx);
        YColNames = [YColNames, 'CV1LCO'];
    end

    cv2lcoIdx = 0; 
    if isfield(origRAND,'CV2LCO')
        cv2lcoColIdx = zeros(1,size(Y,2));
        Y = [Y, origRAND.CV2LCO.ind];
        cv2lcoColIdx =[cv2lcoColIdx, 1];
        cv2lcoIdx = length(cv2lcoColIdx);
        YColNames = [YColNames, 'CV2LCO'];
    end
    

    % check if sample size dependent vectors are defined in preproc 
    % - estimate betas in subgroup only
    % - ... ???
    preprocVecs.ModPStep = [];
    preprocVecs.Sequence = [];

    preprocIdx = zeros(1,length(Y));
    
    for m = length(mods)
        preprocs = xNM.analysis{1,curanal}.params.TrainParam.PREPROC{1,m}.ACTPARAM;
        for ps = 1:length(preprocs)
            switch preprocs{1,ps}.cmd
                case 'correctnuis' % type of preproc vec param 1 
                    if isfield(preprocs{1,ps},'SUBGROUP')
                        preprocVecs.ModPStep = [preprocVecs.ModPStep; m, ps]; 
                        Y = [Y, preprocs{1,ps}.SUBGROUP];
                        preprocIdx = [preprocIdx, 1];
                        preprocVecs.Sequence = [preprocVecs.Sequence, 1]; 
                        YColNames = [YColNames, "PREPROCPARAM"+m*ps];
                        break;
                    end
                case 'rank' % type 2
                    if isfield(preprocs{1,ps}, 'LABEL')
                        preprocVecs.ModPStep = [preprocVec.ModPStep; m, ps]; 
                        Y = [Y, preprocs{1,ps}.LABEL];
                        preprocIdx = [preprocIdx, 1];
                        preprocVecs.Sequence = [preprocVecs.Sequence, 2];   
                        YColNames = [YColNames, "PREPROCPARAM"+m*ps];
                        break;
                    end
                case 'variance' % TO DO
                    if isfield(preprocs{1,ps}, 'LABEL')
                        preprocVecs.ModPStep = [preprocVec.ModPStep; m, ps];  
                        Y = [Y, preprocs{1,ps}.LABEL];
                        preprocIdx = [preprocIdx, 1];
                        preprocVecs.Sequence = [preprocVecs.Sequence, 2];   
                        YColNames = [YColNames, "PREPROCPARAM"+m*ps];
                        break;
                    end
            end
           
        end
    end

    % ADD RANK, VARIANCE


%    
    Ytab = array2table(Y, 'VariableNames', YColNames);
    writetable(Ytab, origDataFile);

    
    condGroupVec = []
    % if certain condition 
    if condIdx >= 0 
        if condIdx == 0 % label
            condName = 'label';
            condIdxY = size(Ytab,2) + 1;
            condGroupVec = repelem(unique(xNM.label),nc);
        elseif condIdx > 0 % one of the covars
            % check whether it's a categorical variable  TO DO 
            condName = xNM.covnames(condIdx);
            [tf, condIdxY] = ismember(condName, COBREgmv.Properties.VariableNames);
            condGroupVec = repelem(unique(xNM.covars(condIdx)),nc);
        end
    end

    
    R = zeros(1,reps); % repeated for 10 times to increase stability of results
    for k=1:reps
        tempRAND = origRAND;
        tic
        M_file = pyrunfile('py_simulate_data.py', 'out_path', ...
            data_file = origDataFile, ...
            labels = double(labels), ...
            n_obs = int64(nc), ... % if vector, then n observations to be simulated within each group (defined by label)
            cv1lco = int64(cv1lcoIdx), ... % vector of 
            cv2lco = int64(cv2lcoIdx), ...
            sitesCols = int64(sitesIdxY), ... % must have length > 1 (dummy coded sites), otherwise discarded
            condVals = int64(condGroupVec), ... % vector that has length = n-to-simulate (nc)
            condColName = condName, ... % column name of the variable that provides groups for special sampling % condCol = int64(condIdxY), ...
            rootdir = analrootdir);
        toc
        M = readtable(py2mat(M_file));

        L = table2array(M(:,end));
%         L(L==0) = -1;
        M = table2array(M(:,1:end-1));

        
        % if indicated, data will be added to original data
        if add2orig
            %Mtab = array2table(M);
            %Otab = Y;
            %Mtab.Properties.VariableNames = Otab.Properties.VariableNames;
            M = [Y;M];
            L = [labels;L];
        end
        
        if sum(preprocIdx)>0
            preprocIdx = logical(preprocIdx);
            for pt = 1:size(preprocVecs.Sequence,2)
                m = preprocVecs.ModPStep(pt,1);
                ps = preprocVecs.ModPStep(pt,2);
                if preprocVecs.Sequence(pt) == 1 % covars subgroup
                    xNM.analysis{1,curanal}.params.TrainParam.PREPROC{1,m}.ACTPARAM{1,ps}.SUBGROUP = M(:,size(M,2)-size(preprocVecs.Sequence,2)+pt);
                elseif preprocVecs.Sequence(pt) == 2 % rank target labels
                    xNM.analysis{1,curanal}.params.TrainParam.PREPROC{1,m}.ACTPARAM{1,ps}.LABEL = M(:,size(M,2)-size(preprocVecs.Sequence,2)+pt);
                elseif preprocVecs.Sequence(pt) == 3 % variance matrix
                    % TO DO 
                end
            end
            M = M(:,~preprocIdx);
        end
       
         
         % TO DO: whether group size relation is correct otherwise
         % potentially use "conditions" of sdv package
        if isfield(tempRAND,'CV2LCO')
            % replace old CV2LCO with new one
            cv2lcoColIdx = logical(cv2lcoColIdx);
            tempRAND.CV2LCO.ind = M(:,cv2lcoColIdx == 1);
            M = M(:,~cv2lcoColIdx);
        end

        if isfield(tempRAND,'CV1LCO')
            % replace old CV1LCO with new one
            cv1lcoColIdx = logical(cv1lcoColIdx);
            tempRAND.CV1LCO.ind = M(:,cv1lcoColIdx == 1);
            M = M(:,~cv1lcoColIdx);
        end

        if isfield(xNM, 'covars')
            % add simulated covariates to xNM structure for preproc
            covColIdx = logical(covColIdx);
            xNM.covars = M(:,covColIdx);
            % remove covariate columns
            M = M(:,~covColIdx);
        end

        % separate different modalities into their own containers in xNM
        for i = 1:length(unique(modColIdx))
            xNM.Y{1,i} = M(:,modColIdx == mods(i));
        end

        nmr = size(M,2);
        

      
        xNM.analysis{1,curanal}.params.TrainParam.FUSION.M = mods;
        xNM.label = L;
%         origCV2LCO = [];

        % what is CV2LCO again?? I assume it does not work if simulated
        % data is added to only one group
%         if add2orig
%             if isfield(RAND,'CV2LCO')
%                 origCV2LCO = sort(RAND.CV2LCO.ind);
%                 groups = unique(origCV2LCO);
%                 nGroups = length(groups);
%                 origN = length(origCV2LCO);
%                 simCV2LCO = [];
%                 for i=1:nGroups
%                     groupN = sum(origCV2LCO == groups(i));
%                     propGroup = groupN/origN;
%                     simGroupN = round(nc*propGroup);
%                     simGroupV = repelem(groups(i), simGroupN);
%                     simCV2LCO = [simCV2LCO, simGroupV];
%                 end
%                 % check whether simCV2LCO == nc, could be slightly different
%                 % due to rounding
%                 if length(simCV2LCO) > nc
%                     % delete one group member for as long as necessary to get
%                     % the vector to the required size
%                     dif = length(simCV2LCO)-nc;
% 
%                     for i = 1:dif
%                         group = groups(i);
%                         groupIDXs = find(simCV2LCO == group);
%                         rmIDX = groupIDXs(end);
%                         simCV2LCO(rmIDX) = [];
%                     end
% 
%                 elseif length(simCV2LCO) < nc
%                     % add one group member for as long as necessary to get
%                     % the vector to the required size
%                     dif = nc-length(simCV2LCO);
%                     for i = 1:dif
%                         group = groups(i);
%                         simCV2LCO(end+1) = group;
%                     end
% 
%                 end

%             end

%             RAND.CV2LCO.ind = simCV2LCO;
%             origCV1LCO = [];
%             if isfield(RAND,'CV1LCO')
%                 origCV1LCO = sort(RAND.CV1LCO.ind);
%                 groups = unique(origCV1LCO);
%                 nGroups = length(groups);
%                 origN = length(origCV1LCO);
%                 simCV1LCO = [];
%                 for i=1:nGroups
%                     groupN = sum(origCV1LCO == groups(i));
%                     propGroup = groupN/origN;
%                     simGroupN = round(nc*propGroup);
%                     simGroupV = repelem(groups(i), simGroupN);
%                     simCV1LCO = [simCV1LCO, simGroupV];
%                 end
%                 % check whether simCV2LCO == nc, could be slightly different
%                 % due to rounding
%                 if length(simCV1LCO) > nc
%                     % delete one group member for as long as necessary to get
%                     % the vector to the required size
%                     dif = length(simCV1LCO)-nc;
% 
%                     for i = 1:dif
%                         group = groups(i);
%                         groupIDXs = find(simCV1LCO == group);
%                         rmIDX = groupIDXs(end);
%                         simCV1LCO(rmIDX) = [];
%                     end
% 
%                 elseif length(simCV1LCO) < nc
%                     % add one group member for as long as necessary to get
%                     % the vector to the required size
%                     dif = nc-length(simCV1LCO);
%                     for i = 1:dif
%                         group = groups(i);
%                         simCV1LCO(end+1) = group;
%                     end
% 
%                 end
% 
%             end
% 
%             RAND.CV1LCO.ind = simCV1LCO;
%         end

        simCV = nk_MakeCrossFolds(L, tempRAND, xNM.modeflag,[], xNM.groupnames, [], 0);
        xNM.cv = simCV;
        xCV = simCV;
        xNM.analind = varargin{1,4};

        xinp = struct('analind',varargin{1,4}, ...
            'lfl',1, ...
            'preprocmat',[], ...
            'gdmat',[], ...
            'gdanalmat', [], ...
            'varstr', [], ...
            'concatfl', [], ...
            'ovrwrt', 2, ...
            'update', true, ...
            'batchflag', true, ...
            'simFlag', true);
        [act, inp] = nk_MLOptimizerPrep(999, [],'Simulate analysis');
        %[params,model] = nk_GetParams2(xNM,)

        %[params,model] = nk_GetParams2(xNM,)
        %R(k) = xNM.analysis{1,varargin{1,4}}.TestPerformanceMean;
        R(k) = xNM.analysis{1,varargin{1,4}}.GDdims{1,1}.BinClass{1,1}.contigency.BAC;
        %         for i = 1:size(xNM.analysis{1,varargin{1,4}}.GDdims{1,1}.predictions{1,1},2)
        %             cur_pred = zeros(size(L));
        %             for pp = 1:size(xNM.analysis{1,varargin{1,4}}.GDdims{1,1}.predictions,1)
        %                 cur_pred(pp,1) = xNM.analysis{1,varargin{1,4}}.GDdims{1,1}.predictions{pp,1}(i);
        %             end
        %             R(k,i) = BAC(L,cur_pred);
        %         end

        fprintf(' ==> mean simulation outcome at %1.3f in repetition: %g,',R(k),k);
    end
else
    nmr 	= ceil(nf*mr);
    nc1     = ceil(nc*er);
    nc2     = nc-nc1;
    verbose = varargin{1,1};
    % Create labels for synthetic data
    L = [ones(nc1,1); -1*ones(nc2,1)];

    IN.ZeroOne = 2;
    FM = [rand(nc1,nmr); -1*rand(nc2,nmr)];
    M = [ FM nk_PerfScaleObj(randn(nc,nf-nmr), IN)];


    % Define algorithm to use for the simulation
    LIBSVMTRAIN = @svmtrain312;
    LIBSVMPREDICT = @svmpredict312;

    switch algorithm
        case 'LINKERNSVM'
            SVM.prog = 'LIBSVM';
            SVM.(SVM.prog).Weighting = 1;
            SVM.(SVM.prog).Optimization.b = 0;
            SCALEMODE = 'scale';
        case 'LINSVM'
            SVM.prog = 'LIBLIN';
            SVM.(SVM.prog).b = 0;
            SVM.(SVM.prog).classifier = 3;
            SVM.(SVM.prog).tolerance = 0.01;
            CMDSTR.WeightFact = 1;
            SCALEMODE = 'std';
        case 'L2LR'
            SVM.prog = 'LIBLIN';
            SVM.(SVM.prog).b = 0;
            SVM.(SVM.prog).classifier = 0;
            SVM.(SVM.prog).tolerance = 0.01;
            CMDSTR.WeightFact = 1;
            SCALEMODE = 'std';
        case 'L1SVC'
            SVM.prog = 'LIBLIN';
            SVM.(SVM.prog).b = 0;
            SVM.(SVM.prog).classifier = 5;
            SVM.(SVM.prog).tolerance = 0.01;
            SCALEMODE = 'std';
        case 'L1LR'
            SVM.prog = 'LIBLIN';
            SVM.(SVM.prog).b = 0;
            SVM.(SVM.prog).classifier = 6;
            SVM.(SVM.prog).tolerance = 0.01;
            SCALEMODE = 'std';
        case 'RF' % not implemented yet
    end

    % We basically always want to weight the decision boundary
    SVM.(SVM.prog).Weighting = 1;
    CMDSTR.WeightFact = 1;

    % At the moment only binary classification is supported
    MODEFL = 'classification';

    % We want to repeat this ten times for producing more stable results
    reps = 10;
    pCV = RAND.OuterPerm;
    nCV = RAND.OuterFold;
    R = zeros(reps,pCV*nCV);

    featmode = 'exponential';
    xSVM = SVM;
    % It is wonderful that we can use multi-core processors nowerdays,
    % (if the Parallel Computation toolbox is available!)
    %%
    xfromData = fromData;
    parfor k=1:reps

        Ystar = M;
        auc_ystar=zeros(1,nmr);

        fprintf(['\n\nRepetition %g:\nWorking on new feature matrix:' ...
            '\n\t%g\tcases' ...
            '\n\t%g\tfeatures' ...
            '\n\t%1.2f\tevent rate', ...
            '\n\t%g\tpredictive markers' ...
            '\n\t%1.2f\tAUCmax' ...
            '\n\t%1.2f\tAUCmin' ...
            '\n\t%g\tbatches', ...
            '\n\t%1.2f\tbatch effect-to-feature scale ratio', ...
            '\n\t%1.2f\tratio of cases with missings', ...
            '\n\t%1.2f\tratio of features with missing per case\n'], ...
            k, nc, nf, er, nmr, auc_max, auc_min, nb, bp, ncm, nfm);

        % CREATE ARTIFICIAL DATA
        if ~xfromData
            switch featmode
                case 'simple'
                    for u=1:nmr
                        fprintf('.')
                        Ystar(:,u) = nk_CreatPredFeat(auc_max, auc_min, L, false);
                    end
                case 'exponential'
                    gran = 5;
                    sepvec = nk_ExpVec(auc_min,auc_max, gran);
                    featvec = floor(linspace(1,nmr,gran));
                    sepvec2 = zeros(nmr,2);
                    for u=1:gran-1
                        sepvec2(featvec(u):featvec(u+1),:) = repmat([sepvec(gran-u) sepvec(gran-u+1)], featvec(u+1)-featvec(u)+1,1);
                    end
                    sepvec2(end,:) = [sepvec(gran-u) sepvec(gran-u+1)];
                    if verbose, fprintf('\nRunning exponential decay function.\n'); end
                    for u=1:nmr
                        fprintf('.')
                        Ystar(:,u) = nk_CreatPredFeat(sepvec2(u,2), sepvec2(u,1), L, false);
                    end
            end


            % ADD BATCH EFFECTS
            if nb
                if verbose, fprintf('\nAdding batch effects to matrix.'); end
                % Create exponential batch vector and offsets
                batchvec = floor(nk_ExpVec(1, nc, nb));
                batchoff = rescale(rand(nb,1), -bp, bp);
                % Create batch matrix
                Ybatch = zeros(nc, nf);
                for u=1:nb-1
                    Ybatch(batchvec(u):batchvec(u+1),:) = rand(batchvec(u+1)-batchvec(u)+1,nf)*batchoff(u);
                end
                Ystar = nk_PerfScaleObj(Ystar);
                % Make sure that the order of subjects is permuted before batch
                % effects are added to the matrix
                Iperm = randperm(nc);
                % add batch effects
                Ystar = Ystar + Ybatch(Iperm,:);
                Lk = L;
            else
                Lk = L;
            end

            % ADD MISSINGS

            if ncm
                if verbose, fprintf('\nInserting missings to matrix.'); end
                % creating boolean vector
                b_cm = false(nc,1);
                Ib_cm = randperm(nc,ceil(ncm*nc));
                b_cm(Ib_cm) = true;
                nbcm = numel(Ib_cm);
                Ymiss = zeros(nbcm,nf);
                for u=1:nbcm
                    b_fm = false(1,nf);
                    Ib_fm = randperm(nf,ceil(nfm*nf));
                    b_fm(Ib_fm) = true;
                    Ymiss(u,b_fm) = NaN;
                end
                Ystar(b_cm,:) = Ystar(b_cm,:) + Ymiss;
            end
        end
        Lk = L;
        if verbose, fprintf('\nFeature matrix creation completed.\n'); end
        for u=1:size(Ystar,2), auc_ystar(u) = AUC(Lk,Ystar(:,u)); end

        %Create CV structure using NM
        %Produce label that NM understands
        Lx = Lk; Lx(Lk==-1)=2;
        cv = nk_MakeCrossFolds(Lx, RAND, 'classification', [], {'A','B'} );
        Ri = zeros(pCV,nCV);

        for j=1:pCV

            % Use multi-core setting of my laptop to parallelize the CV cycle
            for i=1:nCV

                Tr = cv.TrainInd{j,i};
                Ts = cv.TestInd{j,i};
                Lr = Lk(Tr);
                Ls = Lk(Ts);
                W = ones(numel(Lr),1);

                % Scale / standardize matrix using NM and do it properly within the cross-validation cycle
                switch SCALEMODE
                    case 'scale'
                        if verbose, fprintf('\nScale data mtrix.'); end
                        [Tr_Ystar, INi] = nk_PerfScaleObj(Ystar(Tr,:)); % Train scaling model
                        Ts_Ystar = nk_PerfScaleObj(Ystar(Ts,:), INi); % Apply it to the test data
                    case 'std'
                        if verbose, fprintf('\nStandardize data mtrix.'); end
                        [Tr_Ystar, INi] = nk_PerfStandardizeObj(Ystar(Tr,:)); % Train standardization model
                        Ts_Ystar = nk_PerfStandardizeObj(Ystar(Ts,:), INi); % Apply it to the test data
                end

                if ncm
                    if verbose, fprintf('\nImpute data matrix.'); end
                    [Tr_Ystar, INi] = nk_PerfImputeObj(Tr_Ystar); % Train imputation model
                    Ts_Ystar = nk_PerfImputeObj(Ts_Ystar, INi); % Apply it to the test data
                end

                switch algorithm
                    case 'LINKERNSVM'
                        % standard LIBSVM params:
                        cmd = '-s 0 -t 0 -c 1';
                        % set weighting!
                        cmd = nk_SetWeightStr(xSVM, MODEFL, CMDSTR, Lr, cmd);
                        % Train model
                        model = LIBSVMTRAIN(W, Lr, Tr_Ystar, cmd);
                        % Test model
                        [ ~, ~, ds ] = LIBSVMPREDICT( Ls, Ts_Ystar, model, sprintf(' -b %g',xSVM.LIBSVM.Optimization.b));
                    case {'LINSVM', 'L2LR', 'L1LR', 'L1SVC'}
                        %Define command string
                        cmd = nk_DefineCmdStr(xSVM, MODEFL);
                        cmd = [ ' -c 1' cmd.simplemodel cmd.quiet ];
                        % set weighting!
                        cmd = nk_SetWeightStr(xSVM, MODEFL, CMDSTR, Lr, cmd);
                        % Train model
                        model = train_liblin244(Lr, sparse(Tr_Ystar), cmd);
                        % Test model
                        [ ~, ~, ds ] = predict_liblin244(Ls, sparse(Ts_Ystar), model, sprintf(' -b %g -q',xSVM.LIBLIN.b));
                        ds = ds(:,1);
                    case 'RF'
                end

                % We use Balanced Accuracy for now.
                Ri(j,i) = BAC(Ls,ds);
            end
        end
        R(k,:) = Ri(:); fprintf(' ==> mean simulation outcome at %1.3f in repetition: %g,',mean(Ri(:)), k);
    end
end

R=R(:);
Rmean = nm_nanmean(R);
R95CI = nm_95confint(R);


