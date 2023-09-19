% =========================================================================
% =                        CLASSIFICATION PLOTS                           =
% =========================================================================
function handles = display_regrplot_oocv_labels_unknown(h, handles)

% Prepare axes
warning off
axes(handles.axes1); cla(handles.axes1); hold on
handles.axes1.Position = handles.axes1pos_orig;
handles.axes38.Visible = 'off';  cla(handles.axes38);
handles.tglSort.Enable =  "off";
MS = 15;
MSoocv = 40;

% Get preferred figure type
GraphType = get(handles.selYaxis,'Value');

[handles, oocvind] = get_oocvind(handles);

% Check whether the labels are known
labels_known = handles.OOCVinfo.Analyses{handles.curranal}.labels_known(oocvind);

clrswp = handles.tglClrSwp.Value;
if clrswp
    g1 = 2; g2 = 1;
else
    g1 = 1; g2 = 2;
end

% Current label
l=1;

% Determine if we have a multi-class analysis

if handles.tglPercRank.Value
    switch GraphType
        case {1,2,3}
            refsample =  handles.Regr.mean_predictions;
            targetsample = handles.OOCV(oocvind).data.RegrResults{l}.MeanCV2PredictedValues;
            zeroline = 0;
        case 4
            refsample =  handles.Regr.prob_predictions(:,1);
            targetsample = handles.OOCV(oocvind).data.RegrResults{l}.BinMajVoteProbabilities;
            zeroline = 0;
        case {5,6}
            refsample =  handles.Regr.mean_predictions;
            targetsample = handles.OOCV(oocvind).data.RegrResults{l}.BinMajVoteProbabilities;
            zeroline = 0;
    end
    zeroln = nk_ComputePercentiles(refsample, zeroline,'inverse');
    P_h = ranktransform(refsample);
    P_oocv_h    = ranktransform(refsample, targetsample);
else
    switch GraphType
    case {1,2,3}
        P_h         = handles.Regr.mean_predictions;
        P_oocv_h    = handles.OOCV(oocvind).data.RegrResults{l}.MeanCV2PredictedValues;
        zeroln = 0;

    end
end

% Get subindex if available
if isfield(handles,'SubIndex'), subfl = true; SubI = handles.SubIndex; else, subfl = false; SubI = true(size(P_oocv_h,1),1); end



legvecn = false(1,2);
handlevecn = cell(1,2);

% Print CV data: Group 1
violin(P_h,'x',1,'facecolor',handles.colptin,'edgecolor','none', 'facealpha', 0.1, 'bw',0.3);% 'mc','','medc','');
handlevecn{1} = dotdensity( 1, P_h,...
    'dotEdgeColor', handles.colptin, ...
    'dotFaceColor',handles.colptin, ...
    'dotSize',MS, ...
    'dotMarker','o', ...
    'dotAlpha', 0.5, ...
    'meanLine', 'off', ...
    'medianLine', 'off'); 

legvecn(1)=1;



if labels_known
   
else
    violin(P_oocv_h,'x',2,'facecolor','k','edgecolor','none', 'facealpha', 0.3, 'bw',0.3);% 'mc','','medc','');
    [handlevecn{2},~,N,X,Y] = dotdensity( 2 , P_oocv_h, ...
        'dotEdgeColor', 'k', ...
        'dotFaceColor', 'k', ...
        'dotSize',MSoocv, ...
        'dotMarker', 'o', ...
        'meanLine', 'off', ...
        'medianLine', 'off');
    legvecn(2)=1;
end

% Define textbox info data 
pss = cell(1,numel(N));psslen=0;
for i=1:numel(pss)
    expgroupi = 'not labeled';

    pss{i} = sprintf(['Subject ID: %s' ...
                    '\nExpected outcome: %s' ...
                    '\nPredicted outcome: %s'], handles.OOCVinfo.Analyses{handles.curranal}.cases{oocvind}{i}, expgroupi,P_oocv_h(i));
    if size(pss{i},2)> psslen, psslen=size(pss{i},2); pssi = i; end
end
hText = uicontrol('Style','text','String', pss{pssi},'FontSize',11, 'Units','normalized', 'Parent', gcf,'Visible','off'); 
figdata.cases       = handles.OOCVinfo.Analyses{handles.curranal}.cases{oocvind}(N);
figdata.x           = X;
figdata.y           = Y;
figdata.patterntext = pss(N);
figdata.parentui    = handles.pnBinary;
figdata.pnpos       = handles.pnBinary.Position;
figdata.figpos      = handles.figure1.Position;
figdata.hPanel      = uipanel('Units','norm', 'Position',hText.Extent, 'BorderType','etchedout', 'BackgroundColor', [.6 .7 .6], 'Visible','off');
figdata.textHdl     = annotation(figdata.hPanel, 'textbox', 'String','', ...
                            'Interpreter','none', ... 
                            'Color', 'black', ...
                            'BackgroundColor',[.6 .7 .6], ...
                            'Position', [0 0 0.99 0.99], ...
                            'EdgeColor', [.6 .7 .6], ...
                            'LineWidth', 0.1, ...
                            'Margin', 5, ...
                            'FitBoxToText','on', ...
                            'Visible','off');
set(handles.axes1,'UserData',figdata);

% Adjust y scaling depending on prob/rank estimates or decision scores
if GraphType >3 || handles.tglPercRank.Value
    probfx = zeroln;
    P_hx = P_h - probfx;
    handles.Regr.P_h = P_hx;
else
    probfx = 0;
    handles.Regr.P_h = P_h;
end

legvecn = logical(legvecn);
LegendStr = { 'Discovery',...
              'OOCV', ...
              'OOCV: unlabeled/not applicable'};
legendvec = LegendStr(legvecn);
handlevec = handlevecn(legvecn==1);
Hvec =[]; for i = 1:numel(handlevec), Hvec = [Hvec handlevec{i}]; end

% Plot binary class deviding line
xlims = numel(legendvec);
xLimits = get(handles.axes1,'XLim'); xLimitsVec = xLimits(1):xLimits(2); 
handles.axes1.XTick = 0.5:1:xlims+0.5; handles.axes1.YGrid = 'off'; handles.axes1.XGrid = 'on';
zeroline = ones(1,numel(xLimitsVec))*probfx;
plot(handles.axes1,xLimitsVec,zeroline,'k--','LineWidth',handles.ZeroLineWidth)
xlim(handles.axes1, [0.5 xlims+0.5]);
ylim(handles.axes1, 'auto');
handles.axes1.XTickLabel = [];

switch GraphType

    case {1,2,3}
        switch handles.params.TrainParam.SVM.prog
            case {'MikSVM','MKLRVM'}
                algostr = 'RVM probability';
            case 'LIBSVM'
                algostr = 'SVM score';
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
            otherwise
                algostr = [handles.params.TrainParam.SVM.prog 'score'];
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

hx(1) = xlabel('Samples'); 
set(hx(1), ...%'FontSize',handles.AxisLabelSize-2,
    'FontWeight',handles.AxisLabelWeight);
hx(2) = ylabel(['Predicted oucome']); 
set(hx(2), ...%'FontSize',handles.AxisLabelSize-2, ...
    'FontWeight',handles.AxisLabelWeight);
handles.legend_classplot = legend(Hvec, legendvec, 'Location','Best','FontSize', 8,'LineWidth',1);%,'FontSize',handles.LegendFontSize); 
legend('boxon')
flg = 'off'; flg2='off';  

handles.pnRocCmds.Visible = flg2;
handles.pnPieCmds.Visible = flg2;
handles.pnContigCmds.Visible = flg;
handles.cmdMetricExport.Visible = flg;
handles.cmdExportAxes20.Visible = flg;
