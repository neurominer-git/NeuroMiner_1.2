% =========================================================================
% =                             REGRESSION PLOT                           =
% =========================================================================
function handles = display_regrplot(handles, markflag, oocvflag)

if ~exist("oocvflag","var") || isempty(oocvflag)
    oocvflag = false;
end

% Get preferred figure type
GraphType = get(handles.selYaxis,'Value');

if ~exist('markflag','var') || isempty(markflag), markflag=false; end

axes(handles.axes1);
uistack(handles.axes1,'top')

if ~oocvflag
    regrplotstr = '';
    pred    = handles.Regr.mean_predictions;
    errbar  = handles.Regr.std_predictions;
    errbarCI1 = handles.Regr.CI1_predictions;
    errbarCI2 = handles.Regr.CI2_predictions;
    ind     = handles.Regr.index_predictions;
    label   = handles.Regr.labels;
    regrplotcl = 'b';
    cla
    hold on
else
    [handles, oocvind] = get_oocvind(handles);
    hold on
    regrplotstr = '_oocv';
    pred    = handles.OOCV(oocvind).data.RegrResults{handles.curlabel}.MeanCV2PredictedValues;
    errbar  = handles.OOCV(oocvind).data.RegrResults{handles.curlabel}.StdCV2PredictedValues;
    errbarCI1  = handles.OOCV(oocvind).data.RegrResults{handles.curlabel}.CICV2PredictedValues(:,1);
    errbarCI2  = handles.OOCV(oocvind).data.RegrResults{handles.curlabel}.CICV2PredictedValues(:,2);
    ind = true(1,height(pred));
    % Check whether the labels are known
    labels_known = handles.OOCVinfo.Analyses{handles.curranal}.labels_known(oocvind);
    if labels_known
        label = handles.OOCVinfo.Analyses{handles.curranal}.label{oocvind};
    else
        label = [];
    end
    regrplotcl = 'r';
end

if ~isfield(handles.Regr,'Xaxis') || isempty(handles.Regr.Xaxis)
    indnan    = ~isnan(label);
    label     = label(indnan);
    lxL = label;
    lxN = 'Observed targets';
else
    indnan  = ~isnan(handles.Regr.Xaxis);
    label  = label(indnan);
    lxL = handles.Regr.Xaxis(indnan);
    lxN = handles.XaxisName;
end

pred = pred(indnan);
errbar = errbar(indnan);
errbarCI1 = errbarCI1(indnan);
errbarCI2 = errbarCI2(indnan);
ind = ind(indnan);
lx = length(label);

if isfield(handles.Regr,'grouping') && markflag
    ngroups = length(unique(handles.Regr.grouping));
    grouping = handles.Regr.grouping;
    markgroups = true;
else
    ngroups = 1; grouping = ones(lx,1); markgroups = false;
end

% Define textbox info data 
findnan = find(indnan);
pss = cell(1,numel(findnan)); psslen=0;
for i=1:numel(findnan)
      pss{i} = sprintf(['Subject ID [%g]: %s' ...
            '\nObserved Target: %g' ...
            '\nPredicted Target: %g\n'], i, handles.subjects{findnan(i)}, label(i), pred(i));
     if size(pss{i},2)> psslen, psslen=size(pss{i},2); pssi = i; end
end
hText = uicontrol('Style','text','String', pss{pssi},'FontSize',11, 'Units','normalized', 'Parent', gcf,'Visible','off'); 
figdata.x = lxL;
figdata.y = pred;
figdata.patterntext = pss;
figdata.parentui = handles.pnBinary;
figdata.hPanel      = uipanel('Units','norm', 'Position',hText.Extent, 'BorderType','etchedout', 'BackgroundColor', [.6 .7 .6], 'Visible','off');
figdata.textHdl     = annotation(figdata.hPanel, 'textbox', 'String','', ...
                            'Interpreter','none', ... %'VerticalAlign', 'Top', ...
                            'Color', 'black', ...
                            'BackgroundColor',[.6 .7 .6], ...
                            'Position', [0 0 0.99 0.99], ...
                            'EdgeColor', [.6 .7 .6], ...
                            'LineWidth', 0.1, ...
                            'Margin', 5, ...
                            'FitBoxToText','on', ...
                            'Visible','off');
set(handles.axes1,'UserData',figdata);

switch GraphType
    
    case 1
        lgstr{1} = '$\mathbf{Median_{OOT}}$';

    case 2
        lgstr{1} = '$\mathbf{Median_{OOT} \pm 95\%CI}$';
        % Mean predictions with [95%] confidence interval
        L = pred - errbarCI1; U = errbarCI2 - pred;
        handles.(['he' regplotstr])  = errorbar(lxL(ind),pred(ind), L(ind), U(ind),'ko','LineWidth',0.5,'MarkerSize',9);
    
    case 3
        lgstr{1} = '$\mathbf{Median_{OOT} \pm SD}$';
        % Mean predictions with standard deviation
        handles.(['he' regplotstr])  = errorbar(lxL(ind),pred(ind),errbar(ind),'ko','LineWidth',0.5,'MarkerSize',9);
end

handles.(['regrplot' regrplotstr]) = scatter(lxL(ind),pred(ind),'Marker','o','MarkerFaceColor',regrplotcl,'MarkerEdgeColor',rgb('LightBlue'),'MarkerFaceAlpha',0.7, 'MarkerEdgeAlpha',0.2,'SizeData',80);

% Mark points according to "grouping"
xhandles.hg = [];
if markgroups
    for i = 1:ngroups
        indg = (ind == 1 & grouping == i);
        handles.hg(i) = plot(lxL(indg),pred(indg),handles.colpt{i},'LineWidth',1,'MarkerSize',9);
    end
    if isfield(handles.Regr,'grouping') && ngroups>1 
        for i=1:ngroups
            lgstr{i+1} = handles.Regr.groupnames{i};
        end
    end
end
if oocvflag
    xrng = [handles.regrplot.XData lxL'];
    yrng = [handles.regrplot.YData pred'];    
else
    xrng = lxL;
    yrng = pred;
end
xstep = range(xrng)/100;
ystep = range(yrng)/100;
rx = (xstep)*5; ry = (ystep)*5;
xLimitsVec = min(xrng):xstep:max(yrng);
xlim([xLimitsVec(1)-rx xLimitsVec(end)+rx]);
yLimitsVec = min(yrng):xstep:max(yrng);
ylim([yLimitsVec(1)-ry yLimitsVec(end)+ry]);

% Add regression line to plot
xy = min(lxL(ind)):(max(lxL(ind))-min(lxL(ind)))/sum(ind):max(lxL(ind));
try % Statistics toolbox available
    [p,s] = polyfit(lxL(ind),pred(ind),1);
    [yhat,dy] = polyconf(p,xy,s,'predopt','curve');
    handles.(['hline' regrplotstr]) = plot(xy,yhat,regrplotcl,'LineWidth',2);
    handles.(['hline_CI' regrplotstr]) = plotshaded(xy,[yhat+dy; yhat-dy],regrplotcl);
    if ~oocvflag
        lgstr{end+1} = '$\mathbf{\hat{y}_{linear}}$';
        lgstr{end+1} = '$\mathbf{\hat{y}_{95\%CI}}$';
    else
        lgstr{end+1} = '$OOCV\mathbf{\hat{y}_{linear}}$';
        lgstr{end+1} = '$OOCV\mathbf{\hat{y}_{95\%CI}}$';
    end
catch % or not
    handles.(['hline' regrplotstr]) = lsline;
    handles.(['hline_CI' regrplotstr]) = [];
    lgstr{end+1} = '$\mathbf{\hat{y}_{linear}}$';
end
% Prepare legend
switch GraphType
    case 1
        hdlvec = [handles.regrplot handles.hline handles.hline_CI ];
        if oocvflag
            hdlvec = [hdlvec handles.regrplot_oocv handles.hline_oocv handles.hline_CI_oocv ];
        end
    case {2,3}
        hdlvec = [handles.he handles.hg handles.hline handles.hline_CI];
        if oocvflag
            hdlvec = [hdlvec handles.he_oocv handles.hline_oocv handles.hline_CI_oocv ];
        end
end
if ~oocvflag
    handles.legend_classplot = legend(hdlvec, lgstr, 'Location','southeast', 'LineWidth', 1, 'Interpreter','latex');
else
    handles.legend_classplot_oocv = legend(hdlvec, [handles.legend_classplot.String lgstr], 'Location','southeast', 'LineWidth', 1, 'Interpreter','latex');
end

% Set descriptions
xlabel(lxN)
ylabel('Predicted targets')

axes(handles.axes5)
if isfield(handles,'txtPerf'); delete(handles.txtPerf); end
    
% Binarize at median
m = nm_nanmedian(label); set(handles.txtBinarize,'String',m);

