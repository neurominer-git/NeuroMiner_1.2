% =========================================================================
% =                             REGRESSION PLOT                           =
% =========================================================================
function handles = display_regrplot(handles, markflag, oocvflag, oocvprepflag, binarizeflag, alphax)

if ~exist("oocvflag","var") || isempty(oocvflag)
    oocvflag = false;
end

if ~exist("oocvprepflag","var") || isempty(oocvprepflag)
    oocvprepflag = false;
end

if ~exist("binarizeflag","var") || isempty(binarizeflag)
    binarizeflag = true;
end

% Get preferred figure type
GraphType = get(handles.selYaxis,'Value');

if ~exist('markflag','var') || isempty(markflag), markflag=false; end

axes(handles.axes1);
uistack(handles.axes1,'top')
label = [];
if ~oocvflag
    regrplotstr = '';
    lgsufstr = 'OOT'; 
    pred    = handles.Regr.mean_predictions;
    errbar  = handles.Regr.std_predictions;
    errbarCI1 = handles.Regr.CI1_predictions;
    errbarCI2 = handles.Regr.CI2_predictions;
    ind     = handles.Regr.index_predictions;
    regr    = handles.Regr;
    label   = handles.labels;
    subjects = handles.subjects;
    regrplotcl = 'b';
    cla
    hold on
else
    [handles, oocvind] = get_oocvind(handles);
    % Check whether the labels are known
    labels_known = handles.OOCVinfo.Analyses{handles.curranal}.labels_known(oocvind);
    hold on
    regrplotstr = '_oocv';
    lgsufstr = 'OOCV';
    if oocvflag && isfield(handles.OOCV.data.RegrResults{1},'Group') && handles.selSubGroupOOCV.Value>1
        groupidx = handles.selSubGroupOOCV.Value-1;
        subjects = handles.OOCV(oocvind).data.tbl.rownames( handles.OOCV(oocvind).data.RegrResults{handles.curlabel}.Group{groupidx}.Index);
        pred    = handles.OOCV(oocvind).data.RegrResults{handles.curlabel}.Group{groupidx}.MeanCV2PredictedValues;
        errbar  = handles.OOCV(oocvind).data.RegrResults{handles.curlabel}.Group{groupidx}.StdCV2PredictedValues;
        errbarCI1  = handles.OOCV(oocvind).data.RegrResults{handles.curlabel}.Group{groupidx}.CICV2PredictedValues(:,1);
        errbarCI2  = handles.OOCV(oocvind).data.RegrResults{handles.curlabel}.Group{groupidx}.CICV2PredictedValues(:,2);
        if isfield(handles.OOCV(oocvind).data.RegrResults{handles.curlabel}.Group{groupidx},'Regr')
            regr   = handles.OOCV(oocvind).data.RegrResults{handles.curlabel}.Group{groupidx}.Regr;
        else
            regr   = [];
        end
        if labels_known
            label = handles.OOCV(oocvind).data.RegrResults{handles.curlabel}.Group{groupidx}.ObservedValues;
        end
   else
        subjects = handles.OOCV(oocvind).data.tbl.rownames;
        pred    = handles.OOCV(oocvind).data.RegrResults{handles.curlabel}.MeanCV2PredictedValues;
        errbar  = handles.OOCV(oocvind).data.RegrResults{handles.curlabel}.StdCV2PredictedValues;
        errbarCI1  = handles.OOCV(oocvind).data.RegrResults{handles.curlabel}.CICV2PredictedValues(:,1);
        errbarCI2  = handles.OOCV(oocvind).data.RegrResults{handles.curlabel}.CICV2PredictedValues(:,2);
        regr    = handles.OOCV(oocvind).data.RegrResults{handles.curlabel}.Regr;
        if labels_known
            label = handles.OOCVinfo.Analyses{handles.curranal}.label{oocvind};
        end
    end
    ind = true(1,height(pred));   
    regrplotcl = 'r';
end

% Update current Regr structure
handles.curRegr = regr;
handles.curRegr.labels = label;
handles.curRegr.tbl_cont.R2                        = handles.curRegr.R2(handles.curlabel);
handles.curRegr.tbl_cont.r                         = handles.curRegr.r(handles.curlabel);
handles.curRegr.tbl_cont.r_95CI_low                = handles.curRegr.r_95CI_low(handles.curlabel);
handles.curRegr.tbl_cont.r_95CI_up                 = handles.curRegr.r_95CI_up(handles.curlabel);
handles.curRegr.tbl_cont.t                         = handles.curRegr.t(handles.curlabel);
handles.curRegr.tbl_cont.p                         = handles.curRegr.p(handles.curlabel);
handles.curRegr.tbl_cont.MAE                       = handles.curRegr.MAE(handles.curlabel);
handles.curRegr.tbl_cont.NRSMD                     = handles.curRegr.NRSMD(handles.curlabel);
if ~isfield(handles.curRegr.tbl_cont,'array')
    arr                                            = cell2mat(struct2cell( handles.curRegr.tbl_cont));
else
    arr = handles.curRegr.tbl_cont.array;
end
handles.curRegr.tbl_cont.rownames                  = fieldnames(handles.curRegr.tbl_cont);
handles.curRegr.tbl_cont.array                     = arr;

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

switch GraphType
    
    case 1
        lgstr{1} = ['$\mathbf{Median^{' lgsufstr '}}$'];

    case 2
        lgstr{1} = ['$\mathbf{Median^{' lgsufstr '} \pm 95\%CI}$'];
        % Mean predictions with [95%] confidence interval
        L = pred - errbarCI1; U = errbarCI2 - pred;
        handles.(['he' regrplotstr]) = errorbar(lxL(ind),pred(ind), L(ind), U(ind),'ko','LineWidth',0.5,'MarkerSize',9);
    
    case 3
        lgstr{1} = ['$\mathbf{Median^{' lgsufstr '} \pm SD}$'];
        % Mean predictions with standard deviation
        handles.(['he' regrplotstr]) = errorbar(lxL(ind),pred(ind),errbar(ind),'ko','LineWidth',0.5,'MarkerSize',9);
end

handles.(['regrplot' regrplotstr]) = scatter(lxL(ind),pred(ind),'Marker','o','MarkerFaceColor',regrplotcl,'MarkerEdgeColor',rgb('LightBlue'),'MarkerFaceAlpha', alphax, 'MarkerEdgeAlpha',0.2,'SizeData',80);

% Add regression line to plot
xy = min(lxL(ind)):(max(lxL(ind))-min(lxL(ind)))/sum(ind):max(lxL(ind));
try % Statistics toolbox available
    [p,s] = polyfit(lxL(ind),pred(ind),1);
    [yhat,dy] = polyconf(p,xy,s,'predopt','curve');
    handles.(['hline' regrplotstr]) = plot(xy,yhat,regrplotcl,'LineWidth',2);
    handles.(['hline_CI' regrplotstr]) = plotshaded(xy,[yhat+dy; yhat-dy],regrplotcl);
    lgstr{end+1} = ['$\mathbf{\hat{y}_{linear}}^{' lgsufstr '}$'];
    lgstr{end+1} = ['$\mathbf{\hat{y}_{95\%CI}}^{' lgsufstr '}$'];
catch % or not
    handles.(['hline' regrplotstr]) = lsline;
    handles.(['hline_CI' regrplotstr]) = [];
    lgstr{end+1} = ['$\mathbf{\hat{y}_{linear}}^{' lgsufstr '}$'];
end

handles.(['regrplot_lgstr' regrplotstr]) = lgstr;

if oocvprepflag, return; end

% Define textbox info data 
findnan = find(indnan);
pss = cell(1,numel(findnan)); psslen=0;
for i=1:numel(findnan)
      pss{i} = sprintf(['Subject ID [%g]: %s' ...
            '\nObserved Target: %g' ...
            '\nPredicted Target: %g\n'], i, subjects{findnan(i)}, label(i), pred(i));
     if size(pss{i},2)> psslen, psslen=size(pss{i},2); pssi = i; end
end
try
hText = uicontrol('Style','text','String', pss{pssi},'FontSize',11, 'Units','normalized', 'Parent', gcf,'Visible','off'); 
catch
    disp(pssi)
end
if oocvflag
    figdata.cases = handles.OOCVinfo.Analyses{handles.curranal}.cases{oocvind};
else
    figdata.cases = handles.NM.cases;
end
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

% Mark points according to "grouping"
handles.hg = [];
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
xLimitsVec = min(xrng):xstep:max(xrng);
xlim([xLimitsVec(1)-rx xLimitsVec(end)+rx]);
yLimitsVec = min(yrng):xstep:max(yrng);
ylim([yLimitsVec(1)-ry yLimitsVec(end)+ry]);

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
    lgstr = handles.regrplot_lgstr ;
else
    lgstr = [handles.regrplot_lgstr handles.regrplot_lgstr_oocv];    
end
if isfield(handles,'legend_classplot')
    delete(handles.legend_classplot);
end
handles.legend_classplot = legend(hdlvec, lgstr, 'Location','southeast', 'LineWidth', 1, 'Interpreter','latex');
HeightScaleFactor = 1.4;
handles.legend_classplot.Position(4) =  handles.legend_classplot.Position(4) * HeightScaleFactor;

% Set descriptions
xlabel(lxN)
ylabel('Predicted targets')

axes(handles.axes5)
if isfield(handles,'txtPerf'); delete(handles.txtPerf); end
  
if binarizeflag
    % Binarize at median
    m = nm_nanmean(label); set(handles.txtBinarize,'String',m);
    handles.curRegr.b_label = label; 
    handles.curRegr.b_label(label>=m) = 1; 
    handles.curRegr.b_label(label<m) = -1;
    handles.curRegr.b_pred = pred; 
    handles.curRegr.b_pred = handles.curRegr.b_pred - m;
    [handles.curRegr.X, ...
     handles.curRegr.Y, ...
     handles.curRegr.T, ...
     handles.curRegr.AUC] = perfcurve2(handles.curRegr.b_label, handles.curRegr.b_pred, 1);
    handles.curRegr.contigmat = ALLPARAM(handles.curRegr.b_label, handles.curRegr.b_pred);
    handles.curRegr.contigmat.BINARIZATION_THRESHOLD = m;
    binarize_regr(handles);
    display_contigmat(handles);
    %handles.axes20.Visible = "on";
    %handles.axes20.Position(4) = 0.28;
end
