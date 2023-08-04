function handles = display_comparator(handles, act)
% =========================================================================
% function handles = display_comparator(handles, act)
% =========================================================================
% Main function to either perform model comparisons (Quade's test for
% multiple models, Wilcoxon's test for binary comparisons) or allow the
% user to interactively plot (sorted) model performances
% _________________________________________________________________________
% (c) Nikolaos Koutsouleris, 09/2022

if ~exist('act','var') || isempty(act), act = 'stats'; end
if ~isfield(handles,'PerfCompWin')
    WinTag = 'PerfCompWin';
    h = findobj('Tag',WinTag);
    if isempty(h)
       handles = init_fig(handles, WinTag, act);
    else
       handles.PerfTab.Win = h; 
    end
end
h = handles.PerfTab.Win ;
set(0,'CurrentFigure',h)

% -------------------------------------------------------------------------
function handles = init_fig(handles, WinTag, ActionMode)

sz = get(0,'ScreenSize');

cnt=0; analyses=[];
                                    
switch ActionMode
    
    case 'stats'
        
        for i=1:numel(handles.NM.analysis)
            if handles.NM.analysis{i}.status
                cnt = cnt+1;
                %analyses{cnt} = handles.NM.analysis{i}.id;
                analysescnt{cnt} = num2str(cnt);
                data{cnt, 1} = handles.NM.analysis{i}.id;
                data{cnt, 2} = i;
                data{cnt, 3} = '';
                data{cnt, 4} = 'no';
                data{cnt, 5} = '';
                data{cnt, 6} = false;
            end
        end
        
        win_wdth = sz(3)/2; win_hght = sz(4)/2; win_x = sz(3)/2 - win_wdth/2; win_y = sz(4)/2 - win_hght/2;
        handles.PerfTab.Win = uifigure('NumberTitle','off', ...
                'Name','NM Statistical Performance Comparator', ...
                'Tag' ,WinTag, ...
                'MenuBar','none', ...
                'Position', [win_x win_y win_wdth win_hght]);
        gl = uigridlayout(handles.PerfTab.Win, [4 4]);
        gl.RowHeight = {25,'1x',25,25};
        multiflag = false; togmult='off'; 
        if numel(unique(handles.NM.label(~isnan(handles.NM.label))))>2, multiflag = true; togmult='on';end
                                    
        handles.PerfTab.addedit        = uieditfield(gl, ...
                                        'Position', [0.05, 0.87, 0.7, 0.05], ...
                                        'Tooltip', sprintf(['Enter additional predictors here that were obtained outside of your NM analysis.' ...
                                                                '\nMake sure that prediction vectors have the same number of cases\n' ...
                                                                'and predictions match NM outputs.']));
        handles.PerfTab.addedit.Layout.Row = 1;
        handles.PerfTab.addedit.Layout.Column = 3;

        handles.PerfTab.addedit_label  = uilabel(gl,'Text',sprintf('Optionally, enter external predictor outputs [ %g cases, variable # of predictors ]', size(handles.NM.cases,1)));
        handles.PerfTab.addedit_label.Layout.Row = 1;
        handles.PerfTab.addedit_label.Layout.Column = [1 2];

        handles.PerfTab.multiflag      = uibutton(gl, 'state', 'Text','Multi-group', 'Enable', togmult, ...
                                        'Tooltip', 'toggle button to compare multi-group or binary classifiers'); 
        handles.PerfTab.multiflag.Layout.Row = 1;
        handles.PerfTab.multiflag.Layout.Column = 4;

        handles.PerfTab.analysisselect = uitable(gl, 'ColumnName', {'Analyses','Idx','Alias','Reference','Order','Select'}, ... 
                                        'ColumnFormat',{'char','numeric','char',{'yes','no'},'char','logical'},...
                                        'ColumnEditable', [false, false, true, true, true, true],...
                                        'ColumnWidth',{400, 50, 200, 50, 50, 'auto'}, ...
                                        'RowName', analyses,...
                                        'Data', data, ...
                                        'CellEditCallback', {@set_order, handles});
        handles.PerfTab.analysisselect.Layout.Row = 2;
        handles.PerfTab.analysisselect.Layout.Column = [1 4];
                                    
        handles.PerfTab.perfplot_radio = uibuttongroup(gl, ...
                                        'Visible','on', 'Units', 'normalized');
        handles.PerfTab.perfplot_radio.Layout.Row = 3;
        handles.PerfTab.perfplot_radio.Layout.Column = [1 4];
              
        handles.PerfTab.perf_radio1    = uicontrol(handles.PerfTab.perfplot_radio, 'style', 'radiobutton', 'Units', 'normalized',...
                                        'String','Violins of performance distributions', ...
                                        'Position',[0.01 0.05 0.25 1]);
                                    
        handles.PerfTab.perf_radio2    = uicontrol(handles.PerfTab.perfplot_radio, 'style', 'radiobutton', 'Units', 'normalized', ...
                                        'String','Violins of one-vs-all performance deltas', ...
                                        'Position',[0.25 0.05 0.25 1]);
                                    
        handles.PerfTab.perf_radio3    = uicontrol(handles.PerfTab.perfplot_radio, 'style', 'radiobutton', 'Units', 'normalized', ...
                                        'String','Box plots of one-vs-one performance deltas', ...
                                        'Position',[0.525 0.05 0.25 1]);
                                    
        handles.PerfTab.perf_radio4    = uicontrol(handles.PerfTab.perfplot_radio, 'style', 'radiobutton', 'Units', 'normalized', ...
                                        'String','Median (IQR) plot ', ...
                                        'Position',[0.8 0.05 0.20 1]);
                                    
        handles.PerfTab.perfplot_radio.Visible = 'on';
                                    
        handles.PerfTab.fileseltext    = uieditfield(gl);
        handles.PerfTab.fileseltext.Layout.Row = 4;
        handles.PerfTab.fileseltext.Layout.Column = [1 2];
                                    
        handles.PerfTab.fileseldlg     = uibutton(gl, ...
                                        'Text','Save as', ...
                                        'ButtonPushedFcn', {@saveas,handles});
        handles.PerfTab.fileseldlg.Layout.Row = 4;
        handles.PerfTab.fileseldlg.Layout.Column = 3;
        
        handles.PerfTab.tabulate       = uibutton(gl,...
                                        'Text','Compare', ...
                                        'FontWeight', 'bold',...
                                        'BackgroundColor', rgb('lightblue'), ...
                                        'ButtonPushedFcn', {@compare_predictors, handles});
        handles.PerfTab.tabulate.Layout.Row = 4;
        handles.PerfTab.tabulate.Layout.Column = 4;
        guidata(handles.figure1,handles);
        
    case 'visual'
        
        for i=1:numel(handles.NM.analysis)
            if handles.NM.analysis{i}.status
                cnt = cnt+1;
                %analyses{cnt} = handles.NM.analysis{i}.id;
                analysescnt{cnt} = num2str(cnt);
                data{cnt, 1} = handles.NM.analysis{i}.id;
                data{cnt, 2} = i;
                data{cnt, 3} = '';
                data{cnt, 4} = false;
            end
        end
        
        handles.PerfTab.Win = figure('NumberTitle','off','Name','NM Visual Comparator [of OOT performance measures]', 'Tag' ,WinTag,'MenuBar','none');
        rowperfs=[];
        handles.PerfTab.analysisselect = uitable('units', 'normalized', 'position', [0.05, 0.110, 0.9, 0.45], ...
                                        'ColumnName', {'Analyses','Idx', 'Alias','Select'}, ... 
                                        'ColumnFormat',{'char', 'numeric', 'char','logical'},...
                                        'ColumnEditable', [false, false, true, true],...
                                        'ColumnWidth',{300, 50, 'auto', 'auto'}, ...
                                        'RowName', analyses,...
                                        'data', data);
        switch handles.modeflag
            case 'classification'
                perfs1 = {'Balanced Accuracy', ...
                            'Accuracy', ...
                            'Sensitivity', ...
                            'Specificity', ...
                            'Positive Predictive Value', ...
                            'Negative Predictive Value', ...
                            'False Positive Rate', ...
                            'Prognostic Summary Index', ...
                            'Positive Likelihood Ratio', ...
                            'Negative Likelihood Ratio', ...
                            'Diagnostic Odds Ratio', ...
                            'Number Needed to Predict', ...
                            'Number Needed to Diagnose', ...
                            'Youden-Index', ...
                            'Matthews Correlation Coefficient', ...
                            'F-Score', ...
                            'G-Mean'}';
                        
                perfs2 = {'BAC', 'acc', 'sens', 'spec', 'PPV', 'NPV', 'FPR', 'PSI', 'pLR', 'nLR', 'DOR', 'NNP', 'NND', 'Youden', 'MCC', 'Fscore', 'Gmean'}';
                perfs3 = repmat({false},numel(perfs2),1);
                perfs4 = repmat({false},numel(perfs2),1);
                perfs5 = {'Red', 'LawnGreen', 'Blue', 'Gold', 'Orange', 'DarkSalmon', 'Gray', 'BurlyWood', 'DarkGreen', 'Turquoise', 'Fuchsia', 'MediumBlue', 'DarkGray', 'DarkRed', 'LightSkyBlue', 'DarkBlue', 'MediumPurple'}';
                perfs = [perfs1 perfs2 perfs3 perfs4 perfs5];
            case 'regression'
        end
        colors = rgb('getcolors');
        handles.PerfTab.perfselect     = uitable('units', 'normalized', 'position', [0.05, 0.58, 0.9, 0.40], ...
                                        'ColumnName', {'Performance Measures','','Select','Separate Axis', 'Color'}, ... 
                                        'ColumnFormat',{'char','char','logical','logical', colors.name'},...
                                        'ColumnEditable', [true, false, true, true, true],...
                                        'ColumnWidth',{300, 0, 'auto', 'auto', 50}, ...
                                        'RowName', rowperfs,...
                                        'data', perfs);
        handles.PerfTab.visualize    = uicontrol('Style','pushbutton', ...
                                        'units','normalized', ...
                                        'Position',[0.80 0.04 0.15 0.06], ...
                                        'String','Visualize', ...
                                        'FontWeight', 'bold',...
                                        'BackgroundColor', rgb('lightgreen'), ...
                                        'Callback', {@visualize_performances, handles});
        handles.PerfTab.sort         = uicontrol('Style','togglebutton', ...
                                        'units','normalized', ...
                                        'Position',[0.65 0.04 0.15 0.06], ...
                                        'String','Sort', ...
                                        'FontWeight', 'normal',...
                                        'BackgroundColor', rgb('lightgrey'));
        handles.PerfTab.line         = uicontrol('Style','togglebutton', ...
                                        'units','normalized', ...
                                        'Position',[0.50 0.04 0.15 0.06], ...
                                        'String','Line plots', ...
                                        'FontWeight', 'normal',...
                                        'BackgroundColor', rgb('lightgrey'));


end
 
function handles = saveas(src, evt, handles)

if ispc 
    ext = '*.xlsx';
else
    ext = '*.csv';
end

[FileName,PathName] = uiputfile(ext,'Save performance table','CompareTable');
handles.PerfTab.fileseltext.Value = fullfile(PathName, FileName);

function handles = set_order(src, evt, handles)

row = evt.Indices(1); 
switch evt.Indices(2) 
    case 3
        col = 6;
        if ~strcmp(src.Data{row,col},'') 
            src.Data{row, col} = true;
        end
        
    case 4
        col = 4;
        vl = src.Data{row,col};
        src.Data(:,4) = {'no'};
        src.Data{row,col}= vl;
        if strcmp(vl,'yes')
            set(handles.PerfTab.perfplot_radio.Children,'Enable','off');
        else
            set(handles.PerfTab.perfplot_radio.Children,'Enable','on');
        end
            
    case 6
        col = 5;
        if evt.NewData && strcmp(src.Data{row,col},'') 
            if row == 1
                src.Data{row,col} = '1';
            else
            %[m,n] = size(src.Data);
                selected = find(cell2mat(src.Data(1:row-1,5)));
                src.Data{row,col} = num2str(numel(selected)+1);
            end
        else
            src.Data{row,col} = '';
        end

end

function handles = compare_predictors(src, evt, handles)

% Check whether path can be created
pth = fileparts(handles.PerfTab.fileseltext.Value);
if isempty(handles.PerfTab.fileseltext.Value) || isempty(pth), errordlg('Provide a valid output path before tabulating the data.'); return; end

if isfield(handles,'curlabel')
    curlabel = handles.curlabel;
else
    curlabel = 1;
end

if ~isfolder(pth)
    [status, msg] = mkdir(pth);
    if status
        rmdir(pth)
    else
        errordlg(msg)
        return
    end
end

AnalysisSelection = cell2mat(handles.PerfTab.analysisselect.Data(:,6));
OrderSelection = str2num(char(handles.PerfTab.analysisselect.Data(:,5)));
RefSelection = strcmp(handles.PerfTab.analysisselect.Data(:,4),'yes');
RefSelectionSel = RefSelection(AnalysisSelection);

if ~any(AnalysisSelection)
    errordlg('You have to select at least one analysis from the list')
    return
end
if sum(AnalysisSelection) ~= numel(OrderSelection),
    errordlg('You have to specify ordering indices for all selected analyses')
end

AnalysisStrings = handles.PerfTab.analysisselect.Data(:,1);
AnalysisIdx     = handles.PerfTab.analysisselect.Data(:,2);
AnalysisAliasStrings = handles.PerfTab.analysisselect.Data(:,3);
I = strcmp(AnalysisAliasStrings,'');
AnalysisAliasStrings(I) = AnalysisStrings(I);
AnalysisAliasStringsSel = AnalysisAliasStrings(AnalysisSelection);
AnalysisIdxSel = cell2mat(AnalysisIdx(AnalysisSelection));

% Analyse cross-validation structures and optimization criteria
a = zeros(1,sum(AnalysisSelection)); nA=numel(a); 
fd = AnalysisIdxSel;
[~,OrderSelection] = sort(OrderSelection,'ascend'); 
fd = fd(OrderSelection); AnalysisAliasStringsSel=AnalysisAliasStringsSel(OrderSelection);
RefSelectionSel = RefSelectionSel(OrderSelection);
CV = handles.NM.analysis{fd(1)}.params.cv;

for i=1:nA
    a(i) = fd(i);
    if i==1
        [ix, jx] = size(handles.NM.analysis{fd(i)}.params.cv.TrainInd);
        [iy, jy] = size(handles.NM.analysis{fd(i)}.params.cv.cvin{1,1}.TrainInd);
        PARAM = handles.NM.analysis{fd(i)}.params.TrainParam.SVM.GridParam;
        [~,~,PARAMFUN] = nk_GetScaleYAxisLabel(handles.NM.analysis{fd(i)}.params.TrainParam.SVM);
    else
        if size(handles.NM.analysis{fd(i)}.params.cv.TrainInd,1) ~= ix 
            errordlg('The cross-validation structures of the selected analyses have an unequal number of CV2 repetitions');
            return
        elseif size(handles.NM.analysis{fd(i)}.params.cv.TrainInd,2) ~= jx
            errordlg('The cross-validation structures of the selected analyses have an unequal number of CV2 folds');
            return
        elseif PARAM ~= handles.NM.analysis{fd(i)}.params.TrainParam.SVM.GridParam
            errordlg('The optimization criteria must be equal across the selected analyses');
        end
    end
end

if ~isempty(handles.PerfTab.addedit.Value)
    AdditionalPredictors = evalin('base', handles.PerfTab.addedit.Value);
    [nPadd, mPadd] = size(AdditionalPredictors);
    if nPadd ~= size(handles.NM.label,1)
        errordlg('External prediction matrices should contain the same number of cases as the label data in your NM structure');
    end
    nanPadd = sum(isnan(AdditionalPredictors),2)>0;
else
    nPadd = size(handles.NM.label,1); mPadd = 0;
    nanPadd = false(nPadd,1);
end

% Determine cases with missing labels or data across selected analyses to
% find a population shared by all predictors
nanL    = sum(isnan(handles.NM.label),2)>0;
nanAnal = false(nPadd, nA);

for i=1:nA
    M = handles.NM.analysis{a(i)}.params.TrainParam.FUSION.M;
    tNanAnal = false(nPadd, numel(M));
    for j=1:numel(M)
        tNanAnal(:,j) = sum(isnan(handles.NM.Y{M(j)}),2) == size(handles.NM.Y{M(j)},2);
    end
    nanAnal(:,i) = any(tNanAnal,2);
end

nanO = any([nanPadd nanL nanAnal],2);
Lg = handles.NM.label(:,curlabel);
Lg(nanO)=NaN;

% Do we have to recompute the prediction performance because of cases with
% NaN labels or complete NaN data? 
if any(nanO), recomp = true; else, recomp = false; end   

nA = sum(AnalysisSelection);
if mPadd>0
    PNames = cellstr([repmat('ExtPred_',mPadd,1) num2str((1:mPadd)')])';
else
    PNames = [];
end

switch handles.PerfTab.multiflag.Value

    case 1 % MULTI-CLASS ANALYSIS
        G    = zeros(ix*jx, mPadd+nA, handles.ngroups);
        G_Multi = zeros(ix*jx, mPadd+nA);
        %Create one-vs-rest labels
        Lgfd = zeros(size(Lg,1),handles.ngroups);
        for curclass=1:handles.ngroups
             ind1 = Lg == curclass;
             ind2 = ~ind1;
             Lgfd(ind1,curclass) = 1; Lgfd(ind2,curclass)=-1;
        end
        Lgfd(nanO,:)=NaN;
        Gnames = cell(ix*jx,handles.ngroups);
        Gnames_Multi = cell(ix*jx,1);
        AnalNames = cell(handles.ngroups,1);
        for curclass=1:handles.ngroups
            
             if mPadd>0
                %Compute performances for external predictors
                for g=1:mPadd
                    ll=1;
                    for f=1:ix
                        for d=1:jx
                             TsInd = CV.TestInd{f,d};
                             G(ll,g,curclass) = PARAMFUN(Lgfd(TsInd), AdditionalPredictors(TsInd));
                             ll=ll+1;
                        end
                    end
                end
             end
             
             lx = size(handles.NM.label,1); ig = mPadd+1;
             % now either get CV2 grids straightaway or recompute grid 
             [ylm, Crit] = nk_GetScaleYAxisLabel(handles.NM.analysis{a(1)}.params.TrainParam.SVM);
             for i=1:nA
                 
                AggrFlag = handles.NM.analysis{a(i)}.params.TrainParam.RFE.CV2Class.EnsembleStrategy.AggregationLevel;
                M = handles.NM.analysis{a(i)}.params.TrainParam.FUSION.M;
                nGDdims = numel(handles.NM.analysis{a(i)}.GDdims);
                
                for g=1:nGDdims

                    AnalG = handles.NM.analysis{a(i)}.GDdims{g};
                    
                    if nGDdims > 1
                        AnalName = sprintf('%s_G%g-vs-REST_M%g', AnalysisAliasStringsSel{i}, curclass, M(g));
                    elseif handles.nclass > 1
                        AnalName = sprintf('%s_G%g-vs-REST', AnalysisAliasStringsSel{i}, curclass);
                    else
                        AnalName = handles.NM.analysis{a(i)}.id;
                    end
                    AnalName = regexprep(AnalName,'-','_');
                    if length(AnalName) > namelengthmax
                        warning('Variable name to long! Removing any underscores');
                        tAnalName = regexprep(AnalName,'_','');
                        if length(tAnalName) > namelengthmax
                            tAnalName = inputdlg(['The variable name is too long (max ' num2str(namelengthmax) ' characters. Please make manual adjustments'],'Error',[], AnalName);
                        end
                        AnalName = tAnalName;
                    end
                    
                    AnalNames{curclass} = [AnalNames{curclass} {AnalName}];

                    if recomp
                        
                        ll=1; 
                        NodesCnt = [ones(lx,1) zeros(lx,1)];

                        for f=1:ix

                            for d=1:jx
                                
                                TsInd = CV.TestInd{f,d};
                                if iscell(handles.NM.analysis{a(i)}.GDdims{g}.multi_bestPpos)
                                    nNodes = numel(handles.NM.analysis{a(i)}.GDdims{g}.multi_bestPpos{ll});
                                else
                                    nNodes = 1;
                                end
                                if ~AggrFlag
                                    nPred = iy*jy;
                                else 
                                    nPred = nNodes*iy*jy;
                                end
                                NodesCnt(TsInd,2) = NodesCnt(TsInd,2) + nPred;
                                llNodesCnt = NodesCnt(TsInd,:); N = numel(TsInd);
                                Pred = AnalG.multi_predictions(TsInd, curlabel); 
                                Pred = arrayfun( @(j) nm_nansum(Pred{j}(llNodesCnt(j,1):llNodesCnt(j,2))==curclass)*100/numel(llNodesCnt(j,1):llNodesCnt(j,2)), 1:N )';
                                indP = Pred > 50;
                                Pred(indP)=1; Pred(~indP)=-1;
                                G(ll,ig,curclass) = PARAMFUN(Lgfd(TsInd,curclass), Pred);
                                if i==1 && g==1
                                    Gnames{ll,curclass} = sprintf('CV2: R%g_F%g_C%g', f,d,curclass); 
                                end
                                ll=ll+1;
                            end
                            NodesCnt(:,1) = NodesCnt(:,2)+1;
                        end
                    else
                        ll=1;
                        for f=1:ix
                            for d=1:jx
                                G(ll,ig,curclass) = AnalG.bestTS{curclass}(f,d);
                                if i==1 && g==1
                                    Gnames{ll,curclass} = sprintf('CV2: R%g_F%g_C%g', f,d,curclass);                                    
                                end
                                ll=ll+1;
                            end
                        end
                    end
                    ig=ig+1;
                end
             end
        end
        % Now work in multi-class mode
        AnalNamesMulti = []; ig = mPadd+1;
        for i=1:nA
            M = handles.NM.analysis{a(i)}.params.TrainParam.FUSION.M;
            nGDdims = numel(handles.NM.analysis{a(i)}.GDdims);
          
            for g=1:nGDdims
                AnalG = handles.NM.analysis{a(i)}.GDdims{g};
                if nGDdims > 1
                    AnalName = sprintf('%s_MultiClass_M%g', AnalysisAliasStringsSel{i}, M(g));
                elseif handles.nclass > 1
                    AnalName = sprintf('%s_MultiClass', AnalysisAliasStringsSel{i});
                else
                    AnalName = handles.NM.analysis{a(i)}.id;
                end
                AnalName = regexprep(AnalName,'-','_');
                if length(AnalName) > namelengthmax
                    warning('Variable name to long! Removing any underscores');
                    tAnalName = regexprep(AnalName,'_','');
                    if length(tAnalName) > namelengthmax
                        tAnalName = inputdlg(['The variable name is too long (max ' num2str(namelengthmax) ' characters. Please make manual adjustments'],'Error',[], AnalName);
                    end
                    AnalName = tAnalName;
                end
                
                AnalNamesMulti = [AnalNamesMulti {AnalName}];
                
                ll=1;
                for f=1:ix
                    for d=1:jx
                         G_Multi(ll,ig) = AnalG.multi_bestTS(f,d);
                         Gnames_Multi{ll} = sprintf('CV2: R%g_F%g', f,d);      
                         ll=ll+1;
                    end
                end
                ig=ig+1;
            end
        end
        AnalNamesMulti = [PNames AnalNamesMulti];
        [ pth ,nam , ext ] = fileparts(handles.PerfTab.fileseltext.Value);
        
        % Perform stats for each one-vs-rest classifier
        for curclass=1:handles.ngroups
            AnalNames{curclass} = [PNames AnalNames{curclass}];
            Filename = fullfile(pth, [nam sprintf('_G%g-vs-REST', curclass) ext]);
            handles.comparator_stats{curclass}.PredictorNames = AnalNames{curclass};
            handles.comparator_stats{curclass}.PredictorPerformances = G(:,:,curclass);
            if numel(AnalNames{curclass})>2
                if ~any(RefSelectionSel)
                    handles.comparator_stats{curclass} = ...
                        quadetest(G(:,:,curclass), Gnames(:,curclass), AnalNames{curclass}, Filename);
                        display_comparator_prepmatplot(G(:,:,curclass), ylm, handles.PerfTab.perfplot_radio, handles.comparator_stats{curclass}, AnalNames{curclass})
                else
                    [handles.comparator_stats{curclass}, handles.comparator_diffs{curclass}] = ... 
                                                  wilcoxon(G(:,RefSelectionSel, curclass)', G(:,~RefSelectionSel, curclass)', ...
                                                  0.05, Gnames(:,curclass), AnalNames{curclass}, Filename);
                    display_comparatot_preprefplot(handles.comparator_diffs{curclass}, Crit, RefSelectionSel, AnalNames{curclass});
                end
            else
                handles.comparator_stats{curclass} = wilcoxon(G(:,1,curclass), G(:,2,curclass), 0.05);
            end
        end
        % Now perform multigroup test
        Filename = fullfile(pth, [nam '_MultiClass' ext]);
        handles.comparator_stats_multi.PredictorNames = AnalNamesMulti;
        handles.comparator_stats_multi.PredictorPerformances = G_Multi;
        if numel(AnalNamesMulti)>2
            if ~any(RefSelectionSel)
                handles.comparator_stats_multi = ...
                                                  quadetest(G_Multi, Gnames_Multi, AnalNamesMulti, Filename);
                display_comparator_prepmatplot(G_Multi, ylm, handles.PerfTab.perfplot_radio, handles.comparator_stats_multi, AnalNamesMulti)
            else
                [handles.comparator_stats_multi, handles.comparator_diffs_multi] = ... 
                                                  wilcoxon(G_Multi(:,RefSelectionSel)', G_Multi(:,~RefSelectionSel)', 0.05, Gnames_Multi, AnalNamesMulti, Filename);
                display_comparatot_preprefplot(handles.comparator_diffs_multi, Crit, RefSelectionSel, AnalNamesMulti);
            end
        else
            handles.comparator_stats_multi = wilcoxon(G_Multi(:,1), G_Multi(:,2), 0.05);
        end
       
    case 0 % BINARY / REGRESSION ANALYSIS
        % for repeated LSO it may make sense to remove a CV2 fold that
        % contains only subjects from one class ( in these cases NM
        % automatically replaces BAC with either sensitivity or
        % specificity ). The statistical analysis would be run on a mixed
        % performance criterion
        col_skip = [];
        % if LGOflag = true then OOT performances are computed for each
        % leave-group-out index, instead of performing the test on CV2
        % performances
        LGOflag = false;
        % Map external predictor to cross-validation structure 
        % (implement mapping to binary dichotomizers in multi-class case)
        for curclass=1:handles.nclass
            G = zeros(ix*(jx-numel(col_skip)),mPadd+nA);
            if mPadd>0
                for g=1:mPadd
                    ll=1;
                    for f=1:ix
                        for d=1:jx
                            [Lgfd, TsInd] = create_CV2labels(handles.modeflag, Lg, CV, f, d, curclass);
                            G(ll,g) = PARAMFUN(Lgfd, AdditionalPredictors(TsInd));
                            ll=ll+1;
                        end
                    end
                end
            end
            lx = size(handles.NM.label,1); ig = mPadd+1;

            % now either get CV2 grids straightaway or recompute grid
          
            AnalNames = [];
            if LGOflag && ~recomp 
                LSO =  handles.NM.analysis{a(1)}.params.TrainParam.RAND.CV2LCO.ind; 
                nLSO = numel(unique(LSO));
                Gnames = cell(nLSO,1);
                G = zeros(nLSO,mPadd+nA);
            else
                Gnames = cell(ix*(jx-numel(col_skip)),1);
            end
            [ylm, Crit] = nk_GetScaleYAxisLabel(handles.NM.analysis{a(1)}.params.TrainParam.SVM);
            % Loop through analyses
            for i=1:nA

                AggrFlag = handles.NM.analysis{a(i)}.params.TrainParam.RFE.CV2Class.EnsembleStrategy.AggregationLevel;
                M = handles.NM.analysis{a(i)}.params.TrainParam.FUSION.M;
                nGDdims = numel(handles.NM.analysis{a(i)}.GDdims);
                
                % Loop through subanalyses (GDdims)
                % if nGDdims > 1 => late-fusion
                for g=1:nGDdims

                    AnalG = handles.NM.analysis{a(i)}.GDdims{g};
                    
                    if nGDdims > 1
                        AnalName = sprintf('%s_Cl%g_M%g', AnalysisAliasStringsSel{i}, curclass, M(g));
                    elseif handles.nclass > 1
                        AnalName = sprintf('%s_Cl%g', AnalysisAliasStringsSel{i}, curclass);
                    else
                        AnalName = AnalysisAliasStringsSel{i};
                    end
                    AnalName = regexprep(AnalName,'-','_');
                    if length(AnalName) > namelengthmax
                        warning('Variable name to long! Removing any underscores');
                        tAnalName = regexprep(AnalName,'_','');
                        if length(tAnalName) > namelengthmax
                            tAnalName = inputdlg(['The variable name is too long (max ' num2str(namelengthmax) ' characters. Please make manual adjustments'],'Error',[], AnalName);
                        end
                        AnalName = tAnalName;
                    end
                    AnalNames = [AnalNames {AnalName}];

                    if recomp
                        
                        % Check if number of predictions is unequal
                        % and not matched to NodesCnt(:,2). Through
                        % an error to notify the user to complete
                        % the analyses properly before comparing
                        % them
                        if numel(unique(cellfun(@(x) size(x,2), AnalG.predictions(:, curclass, curlabel))))>1
                            error(['\nFound an equal number of predictions across subjects in analysis %g.' ...
                                   '\nCheck whether you have completed the analysis properly.'], a(i));
                        end

                        ll=1; NodesCnt = [ones(lx,1) zeros(lx,1)];
                        for f=1:ix
                            for d=1:jx
                                % Create row header
                                if i==1 && g==1, Gnames{ll} = sprintf('CV2: R%g_F%g', f,d); end
                                % Create labels for CV2 partition
                                [Lgfd, TsInd] = create_CV2labels(handles.modeflag, Lg, CV, f, d, curclass);
                                % Compute predictions for CV2 partition
                                Pred = create_CV2predictions(AnalG, AggrFlag, NodesCnt, TsInd, iy, jy, curclass, curlabel, ll);
                                % Compute performance for CV2 partition
                                G(ll,ig) = PARAMFUN(Lgfd, Pred);
                                ll=ll+1;
                            end
                            NodesCnt(:,1) = NodesCnt(:,2)+1;
                        end
                    else
                        if LGOflag
                            for f = 1:nLSO
                                Perfs = AnalG.BinClass{curclass}.mean_predictions( LSO == f );
                                Labels = Lg ; Labels(Lg == curclass)=1; Labels(Lg ~= curclass) = -1; Labels = Labels( LSO == f );
                                G(f,ig) = feval(Crit, Labels, Perfs);
                                Gnames{f} = sprintf('Group %g', f);
                            end
                        else
                            ll=1;
                            for f=1:ix
                                for d=1:jx
                                    if ~isempty(col_skip) && d==col_skip, continue; end
                                    G(ll,ig) = AnalG.bestTS{curclass}(f,d);
                                    Gnames{ll} = sprintf('CV2: R%g_F%g', f,d);
                                    ll=ll+1;
                                end
                            end
                        end
                    end
                    ig=ig+1;
                end
                if nGDdims>1 && isfield(handles.NM.analysis{a(i)},'META')
                     AnalNames = [AnalNames {sprintf('%s_Cl%g_Bagged', AnalysisAliasStringsSel{i}, curclass)}];
                     ll=1; NodesCnt = [ones(lx,1) zeros(lx,1)];
                     for f=1:ix
                        for d=1:jx
                            % Create labels for CV2 partition
                            [Lgfd, TsInd] = create_CV2labels(handles.modeflag, Lg, CV, f, d, curclass);
                            % Compute predictions for CV2 partition
                            Pred = zeros(size(TsInd,1),nGDdims);
                            for gx = 1:nGDdims
                                AnalG = handles.NM.analysis{a(i)}.GDdims{gx};
                                Pred(:,gx) = create_CV2predictions(AnalG, AggrFlag, NodesCnt, TsInd, iy, jy, curclass, curlabel, ll);
                            end
                            % Compute performance for CV2 partition
                            G(ll,ig+1) = PARAMFUN(Lgfd, nm_nanmean(Pred,2));
                            ll=ll+1;
                        end
                        NodesCnt(:,1) = NodesCnt(:,2)+1;
                    end
                end
            end
            AnalNames = [PNames AnalNames];
            if handles.nclass > 1
                [ pth ,nam , ext ] = fileparts(handles.PerfTab.fileseltext.Value);
                Filename = fullfile(pth, [nam sprintf('_Cl%g ', curclass) ext]);
            else
                Filename = handles.PerfTab.fileseltext.Value;
            end
            handles.comparator_stats{curclass}.PredictorNames = AnalNames;
            handles.comparator_stats{curclass}.PredictorPerformances = G;
            
            if numel(AnalNames)>2
                if ~any(RefSelectionSel)
                    % Run quade test if no reference population has been
                    % defined, thus each model is compared to all other
                    % models
                    handles.comparator_stats{curclass} = quadetest(G, Gnames, AnalNames, Filename);
                    display_comparator_prepmatplot(G, ylm, handles.PerfTab.perfplot_radio, handles.comparator_stats{curclass}, AnalNames)
                else
                    % Compare each model against a reference model
                    AnalNames = [AnalNames(RefSelectionSel) AnalNames(~RefSelectionSel)];
                    [handles.comparator_stats{curclass}, handles.comparator_diffs{curclass}] = wilcoxon(G(:,RefSelectionSel)', G(:,~RefSelectionSel)', 0.05, Gnames, AnalNames, Filename);
                    display_comparatot_preprefplot(handles.comparator_diffs{curclass}, Crit, RefSelectionSel, AnalNames);
                end
            else
                handles.comparator_stats{curclass} = wilcoxon(G(:,1)', G(:,2)', 0.05, Gnames, AnalNames, Filename);
            end

        end
end

function SelectedObjectIndex = DetermineSelectedRadioButton(handles)
cnt=1;
for i=numel(handles.Children):-1:1
    Obj{cnt} = handles.Children(i).String;
    cnt=cnt+1;
end
SelectedObjectIndex = find(strcmp(Obj,handles.SelectedObject.String));

function handles = visualize_performances(src, evt, handles)

lw = 2;
mk = 12;

AnalysisSelection = find(cell2mat(handles.PerfTab.analysisselect.Data(:,4)));
AnalysisStrings = handles.PerfTab.analysisselect.Data(:,1);
AnalysisIdx     = handles.PerfTab.analysisselect.Data(:,2);
AnalysisAliasStrings = handles.PerfTab.analysisselect.Data(:,3);
I = strcmp(AnalysisAliasStrings,'');
AnalysisAliasStrings(I) = AnalysisStrings(I);
AnalysisIdxSel     = cell2mat(AnalysisIdx(AnalysisSelection));

PerfSelection = cell2mat(handles.PerfTab.perfselect.Data(:,3));
PerfFullStrings = handles.PerfTab.perfselect.Data(:,1);
PerfStrings = handles.PerfTab.perfselect.Data(:,2);
PerfSeparate = cell2mat(handles.PerfTab.perfselect.Data(:,4));
PerfColors = handles.PerfTab.perfselect.Data(:,5);

nA = numel(AnalysisSelection);
nP = numel(PerfSelection);
nPs = sum(PerfSelection);
nS = sum(PerfSeparate);

Px = zeros(nA, nPs);
cnt_i=0; 
curclass = 1;
gddims = 1;
fd = AnalysisIdxSel;
for i=1:nA
    cnt_i=cnt_i+1;
    cnt_j=0;
    for j=1:nP
        if PerfSelection(j)
            cnt_j=cnt_j+1;
            if ~isfield(handles.NM.analysis{fd(i)},'GDdims')
                error('Analysis %s not completed!', handles.NM.analysis{fd(i)}.id);
            else
                switch handles.NM.modeflag
                    case 'classification'
                       Px(cnt_i,cnt_j) = handles.NM.analysis{fd(i)}.GDdims{gddims}.BinClass{curclass}.prob_contigency.(PerfStrings{j});
                    case  'regression'
                       Px(cnt_i,cnt_j) = handles.NM.analysis{fd(i)}.GDdims{gddims}.Regr.(PerfStrings{j});
                end
            end
        end
    end
    
end

AnalysisAliasStringsSel = AnalysisAliasStrings(AnalysisSelection);
% Eventually sort data according to means of rows
if handles.PerfTab.Win.Children(2).Value
    if nPs>1
        mPx = mean(Px,2);
    else
        mPx = Px;
    end
    [~,sI] = sort(mPx,'ascend');
    Px = Px(sI,:);
    AnalysisAliasStringsSel = AnalysisAliasStringsSel(sI);    
end

figure; hold on
if ~handles.PerfTab.Win.Children(1).Value
    if nPs>1
        if nS>0
            yyaxis left; 
            h{1} = bar(Px(:,~PerfSeparate(PerfSelection)),'grouped');
            yyaxis right;
            h{2} = bar(Px(:,PerfSeparate(PerfSelection)),'grouped');
        else
            h = bar(Px,'grouped');
        end
    else
        h = bar(Px);
    end
    ax = gca;
else
    if nS>0
        indLC = find(PerfSelection & ~PerfSeparate);
        indL = find(~PerfSeparate(PerfSelection)); 
        indRC = find(PerfSelection & PerfSeparate);
        indR = find(PerfSeparate(PerfSelection)); 
        yyaxis left
        for i=1:numel(indL)
            h{1,i} = plot(1:size(Px,1),Px(:,indL(i)),'-o', ...
                'Color', rgb(PerfColors{indLC(i)}), 'MarkerFaceColor', rgb(PerfColors{indLC(i)}), 'MarkerEdgeColor', 'w', 'LineWidth',lw,'MarkerSize',mk);
        end
        yyaxis right
        for i=1:numel(indR)
            h{2,i} = plot(1:size(Px,1),Px(:,indR(i)),'-s', ...
                'Color', rgb(PerfColors{indRC(i)}),'MarkerFaceColor', rgb(PerfColors{indRC(i)}), 'MarkerEdgeColor', 'w','LineWidth',lw,'MarkerSize',mk);
        end
    else
        indSel = find(PerfSelection);
        for i=1:numel(indSel)
            h{i} = plot(1:size(Px,1),Px(:,i),'-o', ...
                'Color',rgb(PerfColors{indSel(i)}), 'MarkerFaceColor', rgb(PerfColors{indSel(i)}), 'MarkerEdgeColor', 'w', 'LineWidth',lw,'MarkerSize',mk);
        end
    end
    ax = gca;
    ax.XTick=1:size(Px,1);
end

legend(PerfFullStrings(PerfSelection),'Location','best'); 
ax.XTickLabel = AnalysisAliasStringsSel;
ax.XTickLabelRotation = 45;
ax.TickLabelInterpreter='none';
ax.XLabel.String = 'Analyses';
ax.XLabel.FontWeight = 'bold';
ax.YAxis(1).Label.String = 'Performance Measure(s)';
ax.YAxis(1).Label.FontWeight = 'bold';
ax.YAxis(1).FontSize = 12;
if numel(ax.YAxis)>1
    ax.YAxis(2).Label.String = 'Performance Measure(s)';
    ax.YAxis(2).Label.FontWeight = 'bold';
    ax.YAxis(2).FontSize = 12;
end
ax.Box = 'on';
ax.YGrid = 'on';
ax.Color = [0.95 0.95 0.95]; 

% _________________________________________________________________________
function [Lgfd, TsInd] = create_CV2labels(modeflag, Lg, CV, f, d, curclass)

% Create labels for the given CV2 partition
switch modeflag
    case 'classification'
        TsInd = CV.TestInd{f,d}(CV.classnew{f,d}{curclass}.ind);
        Lgfd = zeros(size(Lg,1),1);
        if numel(CV.classnew{f,d}{curclass}.groups)>1
            ind1 = Lg==CV.classnew{f,d}{curclass}.groups(1);
            ind2 = Lg==CV.classnew{f,d}{curclass}.groups(2);
        else
            ind1 = Lg==CV.classnew{f,d}{curclass}.groups;
            ind2 = ~ind1;
        end 
        Lgfd(ind1) = 1; Lgfd(ind2)=-1;
    case 'regression'
        TsInd = CV.TestInd{f,d};
        Lgfd = Lg;
end
Lgfd = Lgfd(TsInd); 

% _________________________________________________________________________
function Pred = create_CV2predictions(AnalG, AggrFlag, NodesCnt, TsInd, iy, jy, curclass, curlabel, ll)

% Evaluate how many hyperparameters nodes have been
% selected at the given CV2 partition
if isfield(AnalG,'bestPpos')
    if iscell(AnalG.bestPpos{curclass})
        nNodes = numel(AnalG.bestPpos{curclass}{ll});
    else
        nNodes = numel(AnalG.bestPpos{curclass}(ll));
    end
else
    nNodes=1;
end

% Aggregate predictions or not?
if ~AggrFlag, nPred = iy*jy; else, nPred = nNodes*iy*jy; end

% Extract predictions and compute median
% prediction for the given CV2 partition
NodesCnt(TsInd,2) = NodesCnt(TsInd,2) + nPred;
llNodesCnt = NodesCnt(TsInd,:); N = numel(TsInd);
PredX = AnalG.predictions(TsInd, curclass, curlabel);
Pred = arrayfun( @(j) nm_nanmedian(PredX{j}(llNodesCnt(j,1):llNodesCnt(j,2))), 1:N )';

% _________________________________________________________________________
function display_comparator_prepmatplot(G, ylm, perfplot_radio, comparator_stats, AnalNames)

mw=[]; sw=[];
rI = DetermineSelectedRadioButton(perfplot_radio);
switch rI
    case {1,4}
        switch rI
            case 1
                D = G;
            case 4
                D = [];  mw = nm_nanmedian(G); sw = abs(mw-percentile(G,25)); sw = [sw;abs(mw-percentile(G,75))];
        end
        str = 'Performance';
        hlinepos = mean(ylm);
    case 2
        D = nk_ComputeMnMdPairWiseDiff(G,'md','meandiff');
        str = 'Mean one-vs.-all \Delta(Performance)';
        hlinepos = 0;
    case 3
        D = nk_ComputeMnMdPairWiseDiff(G,'md','alldiff');
        str = 'One-vs.-one \Delta(Performance)';
        hlinepos = 0;
end
if ~isfield(comparator_stats,'tbl_p_fdr_posthoc')
    warndlg(sprintf('Quade test was not significant (W=%1.2f; P=%0.3f). Skipping visual post-hoc analysis.', comparator_stats.W, comparator_stats.p));
else
    display_classcomparison_matrix(comparator_stats.tbl_p_fdr_posthoc, AnalNames, mw, sw, [], D, hlinepos, str);
end

function display_comparatot_preprefplot(comparator_diffs, Crit, RefSelectionSel, AnalNames)

figure; ax = axes; hold on;
d = comparator_diffs';
mw = nm_nanmedian(d); %sw = abs(mw-percentile(d,5)); sw = [sw;abs(mw-percentile(d,95))];
sw = nm_95confint(d);
bar(ax, mw, 'FaceColor', rgb('LightSteelBlue'));
errorbar(ax, 1:numel(mw), mw, sw(1,:), sw(2,:), 'LineStyle', 'none', 'Color', 'k'); 
ax.XTick = 1:numel(mw);
ax.XTickLabel = AnalNames(~RefSelectionSel);
ax.XTickLabelRotation = 45;
ax.XLim = [0.25 numel(mw)+0.75];
ax.Box='on';
ax.TickLabelInterpreter='none';
ax.FontWeight = 'bold';
ax.FontSize = 14;
ax.YAxis.Label.String = sprintf('Difference: %s',Crit);
