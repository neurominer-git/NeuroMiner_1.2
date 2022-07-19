% =========================================================================
% =                        CLASSIFICATION PLOTS                           =
% =========================================================================
function handles = display_classplot(h, handles)

axes(handles.axes1); cla; hold on

switch handles.tglSort.Value
    case 0
        MS = handles.DataMarkerSize; SrtStr = '';
    case 1
        MS = 6; SrtStr = 'Sorted ';
end

%% Display classification plot
% Define X axis data
if ~isfield(handles.BinClass{h},'Xaxis') || isempty(handles.BinClass{h}.Xaxis), 
    lxL = 1:length(handles.BinClass{h}.labelh);
    lxN = [SrtStr 'Subject No.'];
    AltAx = false;
    set(handles.axes1,'Position', handles.axes1pos_orig)
    set(handles.axes38,'Visible','off'); cla(handles.axes38);
else
    lxL = handles.BinClass{h}.Xaxis;
    lxN = handles.XaxisName;
    AltAx = true;
    set(handles.axes1,'Position', handles.axes1pos_alt);
    set(handles.axes38,'Visible','on')
end

% Get preferred figure type
GraphType = get(handles.selYaxis,'Value');

% Plot errobars if needed
switch GraphType
    case {1,2,3}
       
        predh = handles.BinClass{h}.mean_predictions;
        if handles.tglPercRank.Value
            zeroline = nk_ComputePercentiles(predh, 0,'inverse');
            predh = ranktransform(predh);
        end
        
        switch GraphType

            % Mean predictions with 95%-CI
            case 2    
                errbarCI2 = handles.BinClass{h}.CI2_predictions;
                errbarCI1 = handles.BinClass{h}.CI1_predictions;
                if handles.tglPercRank.Value
                    errbarCI2 = ranktransform(handles.BinClass{h}.mean_predictions, errbarCI2);
                    errbarCI1 = ranktransform(handles.BinClass{h}.mean_predictions, errbarCI1);
                end
                switch handles.tglSort.Value
                    case 0 
                        L = predh - errbarCI1; U = errbarCI2 - predh;
                        handles.classplot = errorbar(lxL, predh, ...
                                            L, U, ...
                                            'k','LineWidth',handles.ErrorMarkerWidth, ...
                                            'MarkerSize',handles.ErrorMarkerSize, ...
                                            'LineStyle','none');
                    case 1
                        handles.classplot = plotshaded(lxL,[errbarCI1'; errbarCI2'], 'k');
                end
            
            % Mean predictions with standard deviation
            case 3
               
                errbar  = handles.BinClass{h}.std_predictions;
                if handles.tglPercRank.Value
                    errbarL = ranktransform(handles.BinClass{h}.mean_predictions, handles.BinClass{h}.mean_predictions-errbar/2);
                    errbarU = ranktransform(handles.BinClass{h}.mean_predictions, handles.BinClass{h}.mean_predictions+errbar/2);
                    errbarL = predh - errbarL;
                    errbarU = errbarU - predh;
                else
                    errbarL = errbar/2;
                    errbarU = errbar/2;
                end
                switch handles.tglSort.Value
                    case 0
                        handles.classplot = errorbar(lxL, predh, ...
                                            errbarL, errbarU, ...
                                            'k','LineWidth',handles.ErrorMarkerWidth, ...
                                            'MarkerSize',handles.ErrorMarkerSize, ...
                                            'LineStyle','none');
                    case 1
                         handles.classplot = plotshaded(lxL,[(predh-errbarL)'; (predh+errbarU)'], 'k');
                end
                       
        end
    
    % Majority voting probabilities
    case 4
        predh = handles.BinClass{h}.prob_predictions(:,1);
        if handles.tglPercRank.Value
            zeroline = nk_ComputePercentiles(predh, 0.5,'inverse');
            predh = ranktransform(predh);
        end

    % Cross-CV2 perm majority voting probabilities (95%-CIs)
    case 5
        predh = handles.BinClass{h}.CV2grid.mean_predictions;
        if handles.tglPercRank.Value
            zeroline = nk_ComputePercentiles(predh, 0.5,'inverse');
            predh = ranktransform(predh);
        end
        errbarCI2 = handles.BinClass{h}.CV2grid.CI2_predictions;
        errbarCI1 = handles.BinClass{h}.CV2grid.CI1_predictions;
        if handles.tglPercRank.Value
            errbarCI2 = ranktransform( handles.BinClass{h}.CV2grid.mean_predictions, errbarCI2);
            errbarCI1 = ranktransform( handles.BinClass{h}.CV2grid.mean_predictions, errbarCI1);
        end
        switch handles.tglSort.Value
            case 0 
                L = predh - errbarCI1; U = errbarCI2 - predh;
                handles.classplot = errorbar(lxL, predh, ...
                                    L, U, ...
                                    'k','LineWidth',handles.ErrorMarkerWidth, ...
                                    'MarkerSize',handles.ErrorMarkerSize, ...
                                    'LineStyle','none');
            case 1
                handles.classplot = plotshaded(lxL,[errbarCI1'; errbarCI2'], 'k');
        end
    
    % Cross-CV2 perm majority voting probabilities with standard deviation
    case 6    
        predh = handles.BinClass{h}.CV2grid.mean_predictions;
        if handles.tglPercRank.Value
            zeroline = nk_ComputePercentiles(predh, 0.5,'inverse');
            predh = ranktransform(predh);
        end
        errbar  = handles.BinClass{h}.CV2grid.std_predictions;
        if handles.tglPercRank.Value
            errbarL = ranktransform(handles.BinClass{h}.mean_predictions, handles.BinClass{h}.CV2grid.mean_predictions-errbar/2);
            errbarU = ranktransform(handles.BinClass{h}.mean_predictions, handles.BinClass{h}.CV2grid.mean_predictions+errbar/2);
            errbarL = predh - errbarL;
            errbarU = errbarU - predh;
        else
            errbarL = errbar/2;
            errbarU = errbar/2;
        end
        switch handles.tglSort.Value
            case 0
                handles.classplot = errorbar(lxL, predh, ...
                                    errbarL, errbarU, ...
                                    'k','LineWidth',handles.ErrorMarkerWidth, ...
                                    'MarkerSize',handles.ErrorMarkerSize, ...
                                    'LineStyle','none');
            case 1
                 handles.classplot = plotshaded(lxL,[(predh-errbarL)'; (predh+errbarU)'], 'k');
        end
end

% Define X axis scaling
if ~AltAx
    r = 1;
    minX = 1; maxX = numel(predh);
    xLimitsVecInfo = (minX:maxX)';
else
    r = (nk_Range(handles.BinClass{h}.Xaxis)/100)*5;
    xLimitsVecInfo = handles.BinClass{h}.Xaxis;
    minX = min(handles.BinClass{h}.Xaxis);
    maxX = max(handles.BinClass{h}.Xaxis);
end
XLIMS = [minX-r maxX+r];
yr = range(predh)*0.05;
YLIMS = [min(predh)-yr max(predh)+yr]; 
xlim(handles.axes1, XLIMS);
ylim(handles.axes1, YLIMS);

% Define textbox info data 
pss = cell(1,numel(predh)); psslen=0;
if handles.BinClass{h}.CoxMode
    offs = handles.BinClass{h}.mean_cutoff_probabilities;
else
    offs = zeros(size(predh,1),1);
end

for i=1:numel(pss)
    if handles.BinClass{h}.labelh(i) > 0
        expgroupi = handles.BinClass{h}.groupnames{1}; 
    else
        expgroupi = handles.BinClass{h}.groupnames{2};
    end
    if predh(i) > offs(i)
       predgroupi = handles.BinClass{h}.groupnames{1}; 
    else
       predgroupi = handles.BinClass{h}.groupnames{2}; 
    end
    if  handles.BinClass{h}.CoxMode
        pss{i} = sprintf(['Subject ID [%g]: %s' ...
                    '\nExpected group: %s' ...
                    '\nPredicted Group: %s' ...
                    '\nScore: %g' ...
                    '\nCutoff: %g'], i, handles.BinClass{h}.cases{i}, expgroupi, predgroupi, predh(i), offs(i));
    else
        pss{i} = sprintf(['Subject ID: %s' ...
                    '\nExpected group: %s' ...
                    '\nPredicted Group: %s' ...
                    '\nScore: %g'], handles.BinClass{h}.cases{i}, expgroupi, predgroupi, predh(i));
    end
    if size(pss{i},2)> psslen, psslen=size(pss{i},2); pssi = i; end
end
hText = uicontrol('Style','text','String', pss{pssi},'FontSize',11, 'Units','normalized', 'Parent', gcf,'Visible','off'); 
figdata.cases       = handles.BinClass{h}.cases;
figdata.x           = xLimitsVecInfo;
figdata.y           = predh;
figdata.patterntext = pss;
figdata.parentui    = handles.pnBinary;
figdata.pnpos       = handles.pnBinary.Position;
figdata.figpos      = handles.figure1.Position;
figdata.hPanel      = uipanel('Units','norm', 'Position',hText.Extent, 'BorderType','etchedout', 'BackgroundColor', [.6 .7 .6], 'Visible','off');
figdata.textHdl     = annotation(figdata.hPanel, 'textbox', 'String','', ...
                            'Interpreter','none', ... %'VerticalAlign', 'Top', ...
                            'Color', 'black', ...
                            'BackgroundColor',[.6 .7 .6], ...
                            'Position', [0 0 0.99 0.99], ...
                            'EdgeColor', [.6 .7 .6], ...
                            'LineWidth', 0.1, ...
                            'Margin', 6, ...
                            'FitBoxToText','on', ...
                            'Visible','off');
                        
set(handles.axes1,'UserData',figdata);

% Mark groups with color
if isfield(handles.BinClass{h},'ind2')
    id1 = handles.BinClass{h}.ind1; 
    id2 = handles.BinClass{h}.ind2; 
else
    id1 = handles.BinClass{h}.ind1; 
    id2 = ~handles.BinClass{h}.ind1; 
end

if size(handles.BinClass{h}.labelh,1) == 1
    labelh = handles.BinClass{h}.labelh'; 
else
    labelh = handles.BinClass{h}.labelh;
end

if GraphType > 3 || handles.tglPercRank.Value
    if handles.tglPercRank.Value
        signpred = sign(predh-zeroline); offy = zeroline;
    else
        signpred = sign(predh-0.5); offy = 0.5;
    end
else
    signpred = sign(predh-offs); offy = 0;
end

err = signpred ~= labelh;
idx1 = id1 & ~err; idx2 = id2 & ~err; b(1) = 0; b(2)= 0;

% Plot Group A results
if sum(idx1)
    if ~AltAx && ~handles.tglSort.Value
        fIdx1 = find(id1);
        v = [ [0 offy]; [0 YLIMS(2)]; [lxL(fIdx1(end))+r YLIMS(2)]; [lxL(fIdx1(end))+r offy] ];             
        patch('Faces',[1 2 3 4], 'Vertices', v, 'FaceColor', handles.colptin(handles.BinClass{h}.groupind(1),:), 'EdgeColor', 'none', 'FaceAlpha', 0.15)
    end
    b(1) = plot(lxL(idx1),predh(idx1),handles.colpt,...
        'MarkerSize',MS,...
        'Color',handles.colptin(handles.BinClass{h}.groupind(1),:),...
        'MarkerFaceColor',handles.colptin(handles.BinClass{h}.groupind(1),:),...
        'LineWidth',handles.DataMarkerWidth,...
        'LineStyle','none');    
else
    b(1) = plot(1,NaN,'LineStyle','none');
end

% Plot Group B results
if sum(idx2)
    if ~handles.BinClass{h}.one_vs_all
        CLP = handles.colpt;
        CLR = handles.colptin(handles.BinClass{h}.groupind(2),:);
    else
        CLP = 'o';
        CLR = rgb('DarkGrey');
    end
    if ~AltAx && ~handles.tglSort.Value
        fIdx2 = find(id2);
        v = [ [lxL(fIdx2(1)) offy]; [lxL(fIdx2(1)) YLIMS(1)]; [XLIMS(2) YLIMS(1)]; [XLIMS(2) offy] ];             
        patch('Faces',[1 2 3 4], 'Vertices', v, 'FaceColor', CLR, 'EdgeColor', 'none', 'FaceAlpha', 0.15)
    end
    b(2) = plot(lxL(idx2),predh(idx2),CLP,...
        'MarkerSize',MS,...
        'Color',CLR,...
        'MarkerFaceColor', CLR,...
        'LineWidth',handles.DataMarkerWidth,...
        'LineStyle','none');
else
    b(2) = plot(1,NaN,'LineStyle','none');
end

% Adjust y scaling depending on prob/rank estimates or decision scores
if GraphType >3 || handles.tglPercRank.Value
    if handles.tglPercRank.Value
        probfx = zeroline;
        predhx = predh-zeroline;
    else
        probfx = 0.5;
        predhx = predh-0.5;
    end
    handles.BinClass{h}.predh = predhx;
else
    probfx = 0;
    handles.BinClass{h}.predh = predh;
end

% Zeroline tweaks
xLimits = get(handles.axes1,'XLim'); xLimitsVec = linspace(xLimits(1),xLimits(2), numel(handles.axes1.XAxis.TickValues)-1);
zeroline = ones(1,numel(xLimitsVec))*probfx;

if handles.BinClass{h}.CoxMode
    zeroline = handles.BinClass{h}.mean_cutoff_probabilities;
    xLimitsVec(1)=[]; xLimitsVec(end)=[];
end

% Error plotting
ide1 = id1 & err; ide2 = id2 & err;
% Mark errors in Group A
x1 = plot(lxL(ide1),predh(ide1), '*', 'Color', handles.colptin(handles.BinClass{h}.groupind(1),:),'MarkerSize',handles.DataMissMarkerSize,'LineWidth',handles.DataMissMarkerWidth);
if handles.BinClass{h}.one_vs_all 
    Color2 = rgb('DarkGrey');
else
    Color2 = handles.colptin(handles.BinClass{h}.groupind(2),:);
end

% Mark errors in Group B
x2 = plot(lxL(ide2),predh(ide2), '*', 'Color', Color2,'MarkerSize',handles.DataMissMarkerSize,'LineWidth',handles.DataMissMarkerWidth);  

% Create legend
switch GraphType
    case {2,3,5,6}
    switch GraphType
        case {2,5}
            errest = '95%-CI';
        case {3,6}
            errest = 'SD';
    end 
        handlevec = [b,x1,x2,handles.classplot];
        legendvec = [handles.BinClass{h}.groupnames(:)',{'misclassified'}, {'misclassified'},{errest}];
    otherwise
        handlevec = [b,x1,x2];
        legendvec = [handles.BinClass{h}.groupnames(:)',{'misclassified'}, {'misclassified'}];
end
handles.axes1.XTickMode='auto'; 
handles.axes1.YGrid='off'; 
handles.axes1.XGrid='off'; 

if AltAx,
    % Display regression lines for alternative X Axis
    xl      = get(gca,'Xlim');
    p       = polyfit(lxL(id1),predh(id1),1);     % p returns 2 coefficients fitting r = a_1 * x + a_2
    [ rho(1), pval(1) ]  = corr(lxL(id1),predh(id1));
    r       = p(1) .* xl + p(2);                   % compute a new vector r that has matching datapoints in x
    hl(1)   = plot(xl,r, '-', 'LineWidth', 2, 'Color', [ handles.colptin(handles.BinClass{h}.groupind(1),:) .5]);
    p       = polyfit(lxL(id2),predh(id2),1);     % p returns 2 coefficients fitting r = a_1 * x + a_2
    [ rho(2), pval(2) ]  = corr(lxL(id2),predh(id2));
    r       = p(1) .* xl + p(2);                   % compute a new vector r that has matching datapoints in x
    hl(2)   = plot(xl,r, '-', 'LineWidth', 2, 'Color', [ handles.colptin(handles.BinClass{h}.groupind(2),:) .5]);
    handlevec = [handlevec hl];
    legendvec = [legendvec, {sprintf('r=%1.2f | p=%.3f', rho(1), pval(1))}, {sprintf('r=%1.2f | p=%.3f', rho(2), pval(2))} ];
    % Display misclassification histogram analysis 
    % Group 1 misclassification histogram
    [err_hist1, Bins1] = hist3([lxL(id1) err(id1)],[10 2]); 
    perr_hist1 = err_hist1(:,2) ./ sum(err_hist1,2);
     % Group 2 misclassification histogram
    err_hist2 = hist3([lxL(id2) err(id2)], Bins1); 
    perr_hist2 = err_hist2(:,2) ./ sum(err_hist2,2);
    axes(handles.axes38); 
    bar(handles.axes38, Bins1{1}, perr_hist1,'BarWidth',1,'FaceColor', handles.colptin(handles.BinClass{h}.groupind(1),:),'FaceAlpha',0.5); 
    xlim(handles.axes38, XLIMS);
    ylim(handles.axes38, [0 1]);
    set(handles.axes38, ...%'FontSize', handles.AxisTickSize, ...
        'FontWeight', handles.AxisTickWeight, ...
        'LineWidth', handles.AxisLineWidth);
    ylabel(handles.axes38,'% Misclassified / Bin');
    hold(handles.axes38,'on'); 
    if handles.BinClass{h}.one_vs_all 
        Color2 = rgb('DarkGrey');
    else
        Color2 = handles.colptin(handles.BinClass{h}.groupind(2),:);
    end
    bar(handles.axes38, Bins1{1}, perr_hist2,'BarWidth',1,'FaceColor', Color2,'FaceAlpha',0.5); 
    hold(handles.axes38,'off'); 
    axes(handles.axes1);
end

plot(xLimitsVec,zeroline,'LineWidth',handles.ZeroLineWidth,'Color',rgb('Grey'))
if handles.BinClass{h}.CoxMode
    %plotshaded(lxL,[(zeroline+(offs/2-zeroline))'; (zeroline-(offs/2-zeroline))'], 'k')
end
switch GraphType

    case {1,2,3}
        switch handles.params.TrainParam.SVM.prog
            case {'MikSVM','MKLRVM'}
                algostr = 'RVM probability';
            case 'LIBSVM'
                if handles.params.TrainParam.SVM.LIBSVM.Optimization.b
                    algostr = 'SVM probability [Platt''s probabilities]';
                elseif isfield(handles.params.TrainParam.SVM,'BBQ') && handles.params.TrainParam.SVM.BBQ.flag==1
                    algostr = 'SVM probability [Bayesian Binning into Quantiles]';
                else
                    algostr = 'SVM decision score';
                end
                
            case 'MVTRVR'
                algostr = 'RVR score';
            case 'MEXELM'
                algostr = 'ELM score';
            case 'LIBLIN'
                switch handles.params.TrainParam.SVM.LIBLIN.classifier
                    case {0,6}
                        algostr = 'LR';
                    otherwise
                        algostr = 'SVM';
                end
                switch handles.params.TrainParam.SVM.LIBLIN.b
                    case 0
                        algostr = [algostr ' Score'];
                    case 1
                        algostr = [algostr ' Probability'];
                end
            case 'matLRN'
                algostr = sprintf('matLearn [ %s ]',char(handles.params.TrainParam.SVM.matLRN.algo));
            case 'WBLCOX'
                algostr = sprintf('Willbur-Cox proportional hazards model: predicted risk');
            otherwise
                algostr = [handles.params.TrainParam.SVM.prog ' score'];
        end
    case 4
        algostr = 'OOT-Probability';
    case 5
        algostr = 'Mean OOT-Probability (95%-CIs)';
    case 6
        algostr = 'Mean OOT-Probability (SD)';
end
if handles.tglPercRank.Value
    algostr = [algostr ' (%-ranks)'];
end
hx(1) = xlabel(lxN); 
set(hx(1), ...%'FontSize',handles.AxisLabelSize-2,
    'FontWeight',handles.AxisLabelWeight);
hx(2) = ylabel([SrtStr algostr]); 
set(hx(2), ...%'FontSize',handles.AxisLabelSize-2, ...
    'FontWeight',handles.AxisLabelWeight);
handles.legend_classplot = legend(handlevec,legendvec, 'Location','Best', 'FontSize', 8,'LineWidth',1);%,'FontSize',handles.LegendFontSize); 
legend('boxon')


%% Display Contigency data
%if isfield(handles,'txtPerf'); delete(handles.txtPerf); end
handles = display_contigmat(handles);

%% Display ROC
if isfield(handles.BinClass{h}.contingency,'AUC')
  flg='on';
    [handles.hroc, handles.hroc_random] = display_roc(handles);
else
  flg = 'off';
end     
handles.axes1.XTickLabelMode = 'auto';
handles.pnRocCmds.Visible = flg;
handles.pnPieCmds.Visible = flg;
handles.pnContigCmds.Visible = flg;
handles.axes20.Visible = flg;
%handles.cmdMetricExport.Visible = flg;
handles.cmdExportAxes20.Visible = flg;

%% Display pie charts
[handles.h1pie, handles.h2pie] = display_piecharts(handles);

%% Display contingency plot
handles.h_contig = display_contigplot(handles, [], handles.BinClass{h}.groupnames);
    
