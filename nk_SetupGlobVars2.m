function [status, paramstr] = nk_SetupGlobVars2(dat, act, dispflag, varind)

global ...
    PREPROC ...
    PARMODE ...
    SVM ...
    RFE ...
    CMDSTR ...
    MULTI ...
    SAV ...
    GRD ...
    MODEFL ...
    SCALE ...
    CV ...
    RVM ...
    RAND ...
    VIS ...
    MKLRVM ...
    DATID ...
    SPM5VER ...
    TRAINFUNC ...
    PREDICTFUNC ...
    EVALFUNC ...
    MLI ...
    OOCV ...
    LIBSVMTRAIN ...
    LIBSVMPREDICT ...
    FUSION ...
    MULTILABEL ...
    VERBOSE ...
    NM ...
    META ...
    TIME ...
    CVPOS ...
    TEMPL ...
    STACKING

paramstr = [];

switch act

    case 'setup_main'

        status = 0;
        try
            if isfield(dat.TrainParam,'LABEL') && ...
                    dat.TrainParam.LABEL.flag && ...
                    strcmp(dat.TrainParam.LABEL.newmode,'classification') && ...
                    dat.TrainParam.RAND.Decompose ~= 9
                CV = adjust_cv(dat);
            else
                CV  = dat.cv;
            end
        catch
            paramstr = 'Cross-validation structure';
        end

        try
            RAND = dat.TrainParam.RAND;
        catch
            paramstr = sprintf('%s\n%s',paramstr,'Cross-validation setup');
        end

        try
            MODEFL  = dat.label.modeflag;
        catch
            paramstr = sprintf('%s\n%s',paramstr,'Type of predictor: Classification / Regression model');
        end

        try
            FUSION   = dat.TrainParam.FUSION;
        catch
            paramstr = sprintf('%s\n%s',paramstr,'FUSION parameters');
        end

        try
            STACKING = dat.TrainParam.STACKING;
        catch
            paramstr = sprintf('%s\n%s',paramstr,'STACKING parameters');
        end

        try
            RAND    = dat.TrainParam.RAND;
        catch
            paramstr = sprintf('%s\n%s',paramstr,'Cross-validation parameters');
        end

        try
            SAV     = dat.TrainParam.SAV;
        catch
            paramstr = sprintf('%s\n%s',paramstr,'Save model settings');
        end

        if isfield(NM,'OOCV')
            try
                OOCV = NM.TrainParam.OOCV;
            catch
                paramstr = sprintf('%s\n%s',paramstr,'Independent test validation parameters');
            end
        end

        try
            MLI = dat.TrainParam.MLI;
        catch
            paramstr = sprintf('%s\n%s',paramstr,'Prediction interpretation parameters');
        end

        if isempty(NM)
            NM = evalin('base','NM');
        end

        if isfield(NM,'time') && ~isempty(NM.time)
            TIME = NM.time;
        end

        dat = NM;

        if size(dat.label,2)>1
            MULTILABEL.flag = true;
            MULTILABEL.dim  = size(dat.label,2);
            if isfield(dat,'labelnames')
                MULTILABEL.desc = dat.labelnames;
            else
                MULTILABEL.desc = cellstr([repmat('L',MULTILABEL.dim,1) num2str((1:MULTILABEL.dim)')]);
            end
        else
            MULTILABEL.flag = false;
            MULTILABEL.dim = 1;
            MULTILABEL.desc = [];
        end

        VERBOSE = dat.TrainParam.verbosity;

        DATID = dat.id;
        SPM5VER = nk_CheckSPMver;

        PARMODE = 0;

        if exist('dispflag','var') && ~isempty(dispflag)
            if ~isempty(paramstr), msgbox(paramstr,'Missing parameters detected!','error'); end
        end

    case 'setup_strat'

        status = 0;
        if exist('varind','var') && ~isempty(varind)

            %             if numel(varind)>1
            %                 error('CRITICAL: Only one modality can be initialized at a time!!!')
            %             elseif isfield(dat,'Y') && varind >= numel(dat.Y);
            %                 error('CRITICAL: Modality index exceeds available modalities !!!')
            %             end

            [PREPROC, ...
                RFE, ...
                GRD, ...
                SVM, ...
                LIBSVMTRAIN, ...
                LIBSVMPREDICT, ...
                RVM, ...
                MKLRVM, ...
                CMDSTR, ...
                MULTI, ...
                VIS, paramstr] = nk_CompatParams2(dat.TrainParam, varind, paramstr);

            if isempty(RFE),        paramstr{end+1} = 'Feature selection parameters'; end
            if isempty(PREPROC),    paramstr{end+1} = 'Preprocessing parameters'; end
            if isempty(GRD),        paramstr{end+1} = 'Grid optimization settings'; end
            if isempty(MULTI),      paramstr{end+1} = 'Multi-group parameters'; end
            if isempty(VIS),        paramstr{end+1} = 'Visualization parameters'; end
            if isfield(dat.TrainParam,'MULTILABEL')
                MULTILABEL.sel = dat.TrainParam.MULTILABEL.sel;
            else
                if MULTILABEL.dim > 1
                    MULTILABEL.sel = true(1,MULTILABEL.dim);
                else
                    MULTILABEL.sel = 1;
                end
            end
            if iscell(PREPROC)
                tPREPROC = PREPROC{1};
            else
                tPREPROC = PREPROC;
            end

            if isfield(tPREPROC,'LABELMOD')
                SCALE.LABELMOD = tPREPROC.LABELMOD;
            else
                SCALE = [];
            end
            [TRAINFUNC, PREDICTFUNC] = nk_DefineTrainPredictFunc(true);
            [~,~,EVALFUNC] = nk_GetScaleYAxisLabel(SVM);

            if isfield(dat.TrainParam,'META')
                META = dat.TrainParam.META;
            else
                META = [];
            end

        end

        if exist('dispflag','var') && ~isempty(dispflag)
            if ~isempty(paramstr), msgbox(paramstr,'Missing parameters detected!','error'); end
        end

    case 'check'

        nvar = numel(dat.Y); paramstr = [];
        status = 0;

        checkfields = {'PREPROC', 'SVM',  'GRD',  'RFE', 'VIS'};
        descriptors = {'Preprocessing', ...
            'ML algorithm', ...
            'ML optimization', ...
            'Feature selection', ...
            'Visualization'};
        if numel(unique(NM.label))>2 && strcmp(NM.modeflag,'classification')
            checkfields = [checkfields 'MULTI'];
            descriptors = [descriptors 'Multi-group'];
        end

        if isfield(dat,'TrainParam') && ~isempty(dat.TrainParam)
            if  nvar > 1,
                if ~isfield(dat.TrainParam,'FUSION') || isempty(dat.TrainParam.FUSION)
                    paramstr = 'Fusion parameters'; status = 1;
                else
                    switch dat.TrainParam.FUSION.flag
                        case {0,1,2}
                            checkfields = [checkfields, 'SAV', 'RAND'];
                            descriptors = [descriptors, 'Saving parameters', 'Cross-validation parameters'];
                            params = dat.TrainParam;
                            [status, paramstr] = check_params(paramstr, params, checkfields, descriptors);
                            if status
                                paramstr = sprintf('Missing Parameters: %s', paramstr(3:end));
                            end
                        case 3
                            paramstr = cell(numel(dat.TrainParam.FUSION.M),1); ll = 1;
                            for i = 1:nvar
                                if ~sum(any(dat.TrainParam.FUSION.M == i)), continue; end
                                params = dat.TrainParam.STRAT{i};
                                [lstatus, paramstr{ll}] = check_params(paramstr{ll}, params, checkfields, descriptors);
                                if lstatus
                                    status = 1; paramstr{ll} = sprintf('Modality #%g: Missing parameters =>%s',i, paramstr{ll}(2:end));
                                else
                                    paramstr{ll} = sprintf('Modality #%g: OK',i);
                                end
                                ll=ll+1;
                            end

                    end

                end
            else
                checkfields = [checkfields, 'SAV', 'RAND'];
                descriptors = [descriptors, 'Saving parameters', 'Cross-validation definitions'];
                params = dat.TrainParam;
                [status, paramstr] = check_params(paramstr, params, checkfields, descriptors);
            end
        else
            paramstr = 'Training parameters'; status = 1;
        end
        if ~status,
            paramstr = [];
        elseif iscell(paramstr),
            paramstr = char(paramstr);
        end
        if ~isfield(dat,'modeflag')
            paramstr = char(paramstr,'Prediction framework undefined'); status = 1;
        end
        if ~isfield(dat,'cv')
            paramstr = char(paramstr,'Cross-validation structure undefined');  status = 1;
        end
        if exist('dispflag','var') && ~isempty(dispflag)
            if ~isempty(paramstr) && status
                msgbox(paramstr,'INCOMPLETE SETUP!','error');
            end
        end

    case 'clear'
        clear   global ...
            PREPROC ...
            PARMODE ...
            SVM ...
            RFE ...
            CMDSTR ...
            MULTI ...
            SAV ...
            GRD ...
            MODEFL ...
            CV ...
            RVM ...
            RAND ...
            VIS ...
            MKLRVM ...
            DATID ...
            SPM5VER ...
            TRAINFUNC ...
            PREDICTFUNC ...
            EVALFUNC ...
            OOCV ...
            MLI ...
            LIBSVMTRAIN ...
            LIBSVMPREDICT ...
            VERBOSE ...
            TEMPL ...
            META ...
            STACKING ...
            CVPOS ...
            TIME
    otherwise
        error(['Option ' act ' not available!'])
end
end

%--------------------------------------------------------------------------
function [status, paramstr] = check_params(paramstr, params, checkfields, descriptors)

istr = []; status = 0;
for i=1:numel(checkfields)
    if ~isfield(params,checkfields{i}) || isempty(params.(checkfields{i}))
        istr = sprintf('%s, %s', istr, descriptors{i}); status = 1;
    end
end
if ~isempty(istr)
    if ~isempty(paramstr)
        paramstr = [paramstr istr];
    else
        paramstr = istr;
    end
end
end

%--------------------------------------------------------------------------
function cvadj = adjust_cv(inp)

cvadj = inp.cv;

ulb = unique(inp.TrainParam.LABEL.newlabel);
if any(~isfinite(ulb))
    NaNflag = true; ind = logical(sum(isfinite(ulb),2));
    ulb = ulb(ind,:);
else
    NaNflag = false;
end
nclass = numel(ulb);

if isfield(inp.TrainParam.LABEL, 'newgroupnames')
    g = inp.TrainParam.LABEL.newgroupnames;
else
    g = [];
end

class_wrapper = @(ind) GenClass(g, ulb, nclass, inp.TrainParam.LABEL.newlabel(ind), inp.TrainParam.RAND.Decompose, NaNflag);
newclass = cellfun(class_wrapper, inp.cv.TrainInd, 'UniformOutput', false);

% class_wrapper2 = @(class) GenClass2(class, )
% newclass1 = cellfun(class_wrapper2, newclass, 'UniformOutput', false);

for i = 1:size(cvadj.class,1)
    for j = 1:size(cvadj.class,2)
        if iscell(cvadj.class{i,j})
            for l = 1:length(cvadj.class{i,j})
                newclass{i,j}{l}.TrainInd = cvadj.class{i,j}{l}.TrainInd; 
                newclass{i,j}{l}.TestInd = cvadj.class{i,j}{l}.TestInd; 
                label_wrapper = @(ind) inp.TrainParam.LABEL.newlabel(ind);
                newclass{i,j}{l}.TrainLabel = cellfun(label_wrapper, newclass{i,j}{l}.TrainInd, 'UniformOutput', false);
                newclass{i,j}{l}.TestLabel = cellfun(label_wrapper, newclass{i,j}{l}.TestInd, 'UniformOutput', false);
            end
%             class_wrapper2 = @(class) GenClass2(class, cvadj.class{i,j}, [], inp.TrainParam.LABEL.newlabel);
%             newclass{i,j} = cellfun(class_wrapper2, newclass, 'UniformOutput', false);
        else
            newclass{i,j}.TrainInd = cvadj.class{i,j}.TrainInd; 
            newclass{i,j}.TestInd = cvadj.class{i,j}.TestInd; 
            class.TrainLabel = inp.TrainParam.LABEL.newlabel(cvadj.class{i,j}.TrainInd);
            class.TestLabel = inp.TrainParam.LABEL.newlabel(cvadj.class{i,j}.TestInd);
%             class_wrapper2 = @(class) GenClass2(class, cvadj.class{i,j}.TrainInd, cvadj.class{i,j}.TestInd, inp.TrainParam.LABEL.newlabel);
%             newclass{i,j} = cellfun(class_wrapper2, newclass, 'UniformOutput', false);
        end
               
    end
end

cvadj.class = newclass;

newclassnew = cellfun(class_wrapper, inp.cv.TestInd, 'UniformOutput', false);

cvadj.classnew = newclassnew;

end

%--------------------------------------------------------------------
function class = GenClass(g, xlb, nclass, lb, decomposeflag, nanflag)

switch decomposeflag

    % Generate binary classification classes (One-Vs-One)
    case 1

        cnt=1; class = cell(nclass*(nclass-1)/2,1);

        for i=1:nclass-1

            ijlb = zeros(1,2);
            ijlb(1) = xlb(i);

            for j=i+1:nclass

                ijlb(2) = xlb(j);
                class{cnt}.groups = ijlb;

                if ~isempty(g{i} ) && ~isempty(g{j})
                    class{cnt}.groupdesc = [g{i} ' vs ' g{j}];
                end

                % Positive label
                ind1 = find( lb == ijlb(1) );
                ind2 = find( lb == ijlb(2) );
                label1 = ones(1,numel(ind1))';
                label2 = -1*ones(1,numel(ind2))';
                class{cnt}.ind = [ind1; ind2];
                class{cnt}.label = [ label1; label2 ];
                if nanflag,
                    indnan = find( ~isfinite(lb) );
                    labelnan = nan(numel(indnan),1);
                    class{cnt}.ind = [class{cnt}.ind; indnan];
                    class{cnt}.label = [ class{cnt}.label; labelnan ];
                end
                cnt=cnt+1;
            end
        end

        % Generate binary classification classes (One-Vs-All)
    case 2
        class = cell(nclass,1);
        for i=1:nclass
            if ~isempty(g{i}), class{i}.groupdesc = [g{i} ' vs ALL']; end
            class{i}.groups(1) = xlb(i);
            indpos = lb==xlb(i); indneg = lb~=xlb(i);
            label = zeros(size(lb));
            if any(indpos),label(indpos) = 1; end
            if any(indneg),label(indneg) = -1; end
            class{i}.label = label;
            class{i}.ind = (1:size(lb,1))';
            if nanflag,
                indnan = find( ~isfinite(lb) );
                labelnan = nan(numel(indnan),1);
                class{i}.ind = [class{i}.ind; indnan];
                class{i}.label = [ class{i}.label; labelnan ];
            end
        end

        % Generate multi-group classification setup
    case 9
        class.groups =  1 : nclass ;
        class.groupdesc = 'Multi-group classification';
        class.label = lb;
        class.ind = (1:size(lb,1))';
        if nanflag,
            indnan = find( ~isfinite(lb) );
            labelnan = nan(numel(indnan),1);
            class.ind = [class.ind; indnan];
            class.label = [ class.label; labelnan ];
        end
end

end

% end
% 
% %--------------------------------------------------------------------------
function class = GenClass2(class, inp, TestInd, Label)
if iscell(class)
    for k = 1:length(class)
        class{k}.TrainInd = inp{k}.TrainInd; 
        class{k}.TestInd = inp{k}.TestInd; 
        label_wrapper = @(ind) Label(ind);
        class{k}.TrainLabel = cellfun(label_wrapper, class{k}.TrainInd, 'UniformOutput', false);
        class{k}.TestLabel = cellfun(label_wrapper, class{k}.TestInd, 'UniformOutput', false);
    end
else
    class.TraiInd = TrainInd; 
    class.TestInd = TestInd;
    class.TrainLabel = Label(TrainInd);
    class.TestLabel = Label(TestInd);
end
end
