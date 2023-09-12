% =========================================================================
% =                        VISUALIZATION PLOT                             =
% =========================================================================
function handles = display_visual(handles)
global st
st.ParentFig    = handles.figure1;
varind          = handles.curmodal;
measind         = handles.selVisMeas.Value;
meas            = handles.selVisMeas.String;
load_selPager(handles)
pageind         = handles.selPager.Value;
page            = handles.selPager.String;
pagemanual      = handles.txtPager.String;
sortfl          = handles.tglSortFeat.Value;
filterfl        = handles.tglVisMeas2.Value;
filterthr       = handles.txtThrVisMeas2.String;

handles.spiderPlotButton.Visible = "on";

axes(handles.axes33); cla; hold on
set(handles.axes33,'TickLabelInterpreter','none')

CorrMatStr = {'Correlation matrix', 'Correlation matrix (P value)','Correlation matrix (P value, FDR)'};
NetworkCorrMatStr = {'Network plot correlation matrix', 'Network plot correlation matrix (P value)', 'Network plot correlation matrix (P value, FDR)'};

v = handles.visdata{varind,handles.curlabel};
if v.params.visflag == 1
    featind = 1:v.params.nfeats;
else
    try
        if ~isempty(pagemanual) 
            featind = eval(pagemanual);
        else
            featind = eval(page{pageind});
        end
    catch
        featind = 1:v.params.nfeats;
    end
end

x = 0.5: numel(featind);
curclass = get(handles.popupmenu1,'Value');

if strcmp(handles.popupmenu1.String{curclass},'Multi-group classifier')
    multiflag = true;
else
    multiflag = false;
end

if measind>numel(meas)
    measind=numel(meas);
    handles.selVisMeas.Value=measind;
end

switch meas{measind} 
    case 'Model P value histogram'
        fl = 'off';
        fl2 = 'off';
    otherwise
        if isfield(v,'PermModel_Crit_Global_Multi') && multiflag
            fl = 'off';
        else
            fl = 'on';
        end
        if any(strcmp(meas{measind}, CorrMatStr))
            fl2 = 'on';
            load_selVisMeas2(handles)
        else
            if any(strcmp(meas{measind}, NetworkCorrMatStr))
                fl2 = 'on';
                load_selVisMeas2(handles)
            else
                fl2 = 'off';
            end
        end
end

switch v.params.visflag
    case 1
        handles.selPager.Enable = "off";
        handles.txtPager.Enable = "off";
        handles.tglSortFeat.Enable = "off";
        handles.cmdExportFeats.Enable = "off";
        handles.selVisMeas2.Enable = "off";
        handles.txtThrVisMeas2.Enable = "off";
        handles.tglVisMeas2.Enable = "off";
        handles.txtAlterFeats.Enable = "off";
        handles.cmdExportFeatFig.Enable = "off";
        handles.spiderPlotButton.Visible = "on";
    otherwise
        handles.selPager.Enable = fl;
        handles.txtPager.Enable = fl;
        handles.tglSortFeat.Enable = fl;
        handles.cmdExportFeats.Enable = fl;
        handles.selVisMeas2.Enable = fl2;
        handles.txtThrVisMeas2.Enable = fl2;
        handles.tglVisMeas2.Enable = fl2;
        handles.txtAlterFeats.Enable = fl;
        handles.cmdExportFeatFig.Enable = fl;
        vlineval = [];
        if strcmp(fl,'on')
            if ~isempty(pagemanual) 
                handles.selPager.Enable = 'off';
            else
                handles.selPager.Enable = 'on';
            end
        end
        handles.spiderPlotButton.Visible = "off";
end

switch meas{measind}
    
    case 'Feature weights [Overall Mean (StErr)]'
        [y, se, miny, maxy, vlineval] = MEAN(v, curclass, multiflag);
        
    case 'Feature weights [Grand Mean (StErr)]'
        [y, se, miny, maxy, vlineval] = MEAN_CV2(v, curclass, multiflag);
            
    case 'CV-ratio of feature weights [Overall Mean]'
        [y, miny, maxy, vlineval] = CVRatio(v, curclass, multiflag);

    case 'CV-ratio of feature weights [Grand Mean]'
        [y, miny, maxy, vlineval] = CVRatio_CV2(v, curclass, multiflag);

    case 'Feature selection probability [Overall Mean]'
        [y, miny, maxy, vlineval] = FeatProb(v, curclass, multiflag);

    case 'Probability of feature reliability (95%-CI) [Grand Mean]'
        [y, miny, maxy, vlineval] = Prob_CV2(v, curclass, multiflag);

    case 'Sign-based consistency'
       [y, miny, maxy, vlineval] = SignBased_CV2(v, curclass, multiflag);

    case 'Sign-based consistency (Z score)'
        [y, miny, maxy, vlineval] = SignBased_CV2_z(v, curclass, multiflag);
        
    case 'Sign-based consistency -log10(P value)'
        [y, miny, maxy, vlineval] = SignBased_CV2_p_uncorr(v, curclass, multiflag);

    case 'Sign-based consistency -log10(P value, FDR)'
        [y, miny, maxy, vlineval] = SignBased_CV2_p_fdr(v, curclass, multiflag);

    case 'Spearman correlation [Grand Mean]'
        [y, miny, maxy, vlineval] = Spearman_CV2(v, curclass, multiflag);

    case 'Pearson correlation [Grand Mean]'
        [y, miny, maxy, vlineval] = Pearson_CV2(v, curclass, multiflag);

    case 'Spearman correlation -log10(P value) [Grand Mean]'
        [y, se, miny, maxy, vlineval] = Spearman_CV2_p_uncorr(v, curclass, multiflag);

    case 'Pearson correlation -log10(P value) [Grand Mean]'
       [y, se, miny, maxy, vlineval] = Pearson_CV2_p_uncorr(v, curclass, multiflag);
       
    case 'Spearman correlation -log10(P value, FDR) [Grand Mean]'
       [y, miny, maxy, vlineval] = Spearman_CV2_p_fdr(v, curclass, multiflag);

    case 'Pearson correlation -log10(P value, FDR) [Grand Mean]'
       [y, miny, maxy, vlineval] = Pearson_CV2_p_fdr(v, curclass, multiflag);

    case 'Permutation-based Z Score [Grand Mean]'
       [y, miny, maxy, vlineval] = PermZ_CV2(v, curclass, multiflag);

    case 'Permutation-based -log10(P value) [Grand Mean]'
       [y, miny, maxy, vlineval] = PermProb_CV2(v, curclass, multiflag);
       
    case 'Permutation-based -log10(P value, FDR) [Grand Mean]'
       [y, miny, maxy, vlineval] = PermProb_CV2_FDR(v, curclass, multiflag);

    case 'Analytical -log10(P Value) for Linear SVM [Grand Mean]'
       [y, miny, maxy, vlineval] = Analytical_p(v, curclass, multiflag);

    case 'Analytical -log10(P Value, FDR) for Linear SVM [Grand Mean]'
        [y, miny, maxy, vlineval] = Analyitcal_p_fdr(v, curclass, multiflag);
        
    case 'Correlation matrix'
        y = v.CorrMat_CV2{curclass};
        miny = min(y); maxy = max(y);

    case 'Correlation matrix (P value)'
        y = v.CorrMat_CV2_p_uncorr{curclass};
        miny = min(y); maxy = max(y);

    case 'Correlation matrix (P value, FDR)'
        y = v.CorrMat_CV2_p_fdr{curclass};
        miny = min(y); maxy = max(y);

    case 'Network plot correlation matrix'
        y = v.CorrMat_CV2{curclass};
        miny = min(y); maxy = max(y);
        %y = graph(y);

    case 'Network plot correlation matrix (P value)'
        y = v.CorrMat_CV2_p_uncorr{curclass};
        miny = min(y); maxy = max(y);
        %y = graph(y);

    case 'Network plot correlation matrix (P value, FDR)'
        y = v.CorrMat_CV2_p_fdr{curclass};
        miny = min(y); maxy = max(y);
        %y = graph(y);

    case 'Model P value histogram'
        if multiflag
            y = nm_nanmean(v.PermModel_Crit_Global); 
            vp = nm_nanmean(v.ObsModel_Eval_Global);
            ve = nm_nanmean(v.PermModel_Eval_Global);
            vp_cv2 = nm_nanmean(v.ObsModel_Crit_CV2(:));
            ve_cv2 = nm_nanmean(v.PermModel_CV2(:));
            vp_ci_cv2 = nm_95confint(v.ObsModel_Crit_CV2(:));
        else
            y = v.PermModel_Crit_Global(curclass,:); 
            vp = v.ObsModel_Eval_Global(curclass);
            ve = v.PermModel_Eval_Global(curclass,:);
            vp_cv2 = nm_nanmean(v.ObsModel_Crit_CV2(:));
            ve_cv2 = nm_nanmean(v.PermModel_CV2(:));
            vp_ci_cv2 = nm_95confint(v.ObsModel_Crit_CV2(:));
        end
        perms = length(v.PermModel_Eval_Global(1,:));
end
        
if multiflag && isfield(v,'PermModel_Crit_Global_Multi')
    if measind == 1
         y = v.PermModel_Crit_Global_Multi; 
         vp = v.ObsModel_Eval_Global_Multi;
         ve = v.PermModel_Eval_Global_Multi;
    else
         y = v.PermModel_Crit_Global_Multi_Bin(measind-1,:); 
         vp = v.ObsModel_Eval_Global_Multi_Bin(measind-1);
         ve = v.PermModel_Eval_Global_Multi_Bin(measind-1,:);
    end
    perms = length(v.PermModel_Eval_Global_Multi);
    meas{measind} = 'Model P value histogram';
end

%if ~(measind == 23 | measind == 24 | measind == 25)
y(~isfinite(y))=0;
%end

switch meas{measind}

    case CorrMatStr
         set(handles.pn3DView,'Visible','off'); set(handles.axes33,'Visible','on'); cla; hold on
         feats = v.params.features; nF = numel(feats);
         if ~isempty(handles.txtAlterFeats.String)
            altfeatsvar = handles.txtAlterFeats.String;
            feats = evalin('base',altfeatsvar);
         end
       
         sI = 1:nF;
         if filterfl
             FltMetr = handles.selVisMeas2.String{handles.selVisMeas2.Value};
             switch FltMetr
                 case 'CV-ratio of feature weights [Overall Mean]'
                     [ythr,~,~,vlineval] = CVRatio(v,curclass,multiflag);
                     if isempty(filterthr)
                         handles.txtThrVisMeas2.String='-2 2'; 
                         filterthr = vlineval;
                     else
                         filterthr = str2double(filterthr);
                     end
                 case 'CV-ratio of feature weights [Grand Mean]'
                     [ythr,~,~,vlineval] = CVRatio_CV2(v,curclass,multiflag);
                     if isempty(filterthr)
                         handles.txtThrVisMeas2.String='-2 2'; 
                         filterthr = vlineval;
                     else
                         filterthr = str2double(filterthr);
                     end
                 case 'Feature selection probability [Overall Mean]'
                     [ythr,~,~,vlineval] = FeatProb(v,curclass,multiflag);
                     if isempty(filterthr)
                         handles.txtThrVisMeas2.String='0.5'; 
                         filterthr = vlineval;
                     else
                         filterthr = str2double(filterthr);
                     end
                 case 'Probability of feature reliability (95%-CI) [Grand Mean]'
                     [ythr,~,~,vlineval] = Prob_CV2(v,curclass,multiflag);
                     if isempty(filterthr)
                         handles.txtThrVisMeas2.String='-0.5 0.5'; 
                         filterthr = vlineval;
                     else
                         filterthr = str2double(filterthr);
                     end
                 case 'Sign-based consistency -log10(P value)'
                     [ythr,~,~,vlineval] = SignBased_CV2_p_uncorr(v,curclass,multiflag);
                     if isempty(filterthr)
                         handles.txtThrVisMeas2.String=num2str(-log10(0.05),'%1.2f'); 
                         filterthr = vlineval;
                     else
                         filterthr = str2double(filterthr);
                     end
                 case 'Sign-based consistency -log10(P value, FDR)'
                     [ythr,~,~,vlineval] = SignBased_CV2_p_fdr(v,curclass,multiflag);
                     if isempty(filterthr)
                         handles.txtThrVisMeas2.String=num2str(-log10(0.05),'%1.2f'); 
                         filterthr = vlineval;
                     else
                         filterthr = str2double(filterthr);
                     end
                 case 'Analytical -log10(P Value) for Linear SVM [Grand Mean]'
                     [ythr,~,~,vlineval] = Analytical_p(v,curclass,multiflag);
                     if isempty(filterthr)
                         handles.txtThrVisMeas2.String=num2str(-log10(0.05),'%1.2f'); 
                         filterthr = vlineval;
                     else
                         filterthr = str2double(filterthr);
                     end
                 case 'Analytical -log10(P Value, FDR) for Linear SVM [Grand Mean]'
                     [ythr,~,~,vlineval] = Analyitcal_p_fdr(v,curclass,multiflag);
                     if isempty(filterthr)
                         handles.txtThrVisMeas2.String=num2str(-log10(0.05),'%1.2f'); 
                         filterthr = vlineval;
                     else
                         filterthr = str2double(filterthr);
                     end
             end
             if sortfl, [~,sI] = sort(abs(ythr),'descend'); end
             y = y(sI,sI);
             feats = feats(sI);
             y = y(featind,featind);
             feats = feats(featind);
             if numel(filterthr)>1
                Ix2 = (ythr >= min(filterthr) | ythr <= max(filterthr))';
             else 
                Ix2 = (abs(ythr)>= filterthr)'; 
             end
             Ix2 = Ix2(sI); Ix2 = Ix2(featind);
             Ix = any(y) & Ix2;
         else
             if sortfl, [~,sI] = sort(mean(y),'descend'); end
             y = y(sI,sI);
             feats = feats(sI);
             y = y(featind,featind);
             feats = feats(featind);
             Ix = any(y);
         end
         y = y(Ix,Ix);
         feats = feats(Ix);
         nF = numel(feats);
         if nF<50, FS = 10; elseif nF<100, FS = 7; elseif nF < 150, FS = 5; else, FS = 0; end
         yy = y; 
         try
             if strcmp(meas{measind}, CorrMatStr{1})
                     yy(itril(size(yy))) = -1.1;
                     imagesc( yy ); cbar = colorbar('Visible','on');
                     clabel = 'Correlation coefficient';
                     handles.axes33.Colormap = customcolormap([0 0.48 0.99 1],[rgb('Red'); rgb('White'); rgb('Blue'); rgb('LightGray'); ]);
                     handles.axes33.CLim = [-1.1 1];
                     cbar.Limits = [-1 1];
             else
                     yy(itril(size(yy))) = 0;
                     imagesc( yy ); cbar = colorbar('Visible','on');
                     if strcmp(meas{measind}, CorrMatStr{3})
                         clabel = 'P value (FDR-corrected)';
                     else
                        clabel = 'P value (uncorrected)';
                     end
                     handles.axes33.Colormap =  customcolormap([0 0.25 0.50 0.75 .99 1],[rgb('Red'); rgb('Yellow'); rgb('Cyan'); rgb('Blue'); rgb('DarkViolet'); rgb('LightGray'); ]);
                     handles.axes33.CLim = [ -log10(0.05) max(yy(:))] ;
                     cbar.Limits = [-log10(0.05) handles.axes33.CLim(end)];
             end
             cbar.Label.String=clabel; cbar.Label.FontWeight='bold';

             if FS>0
                handles.axes33.XAxis.FontSize = FS; 
                handles.axes33.YAxis.FontSize = FS; 
                handles.axes33.XTick = 1:sum(Ix); handles.axes33.XTickLabel = feats; handles.axes33.XTickLabelRotation=45;
                handles.axes33.YTick = 1:sum(Ix); handles.axes33.YTickLabel = feats; 
                handles.axes33.XAxis.Label.FontSize = FS;
                handles.axes33.YAxis.Label.FontSize = FS;
             else
                handles.axes33.XTickMode = 'auto';
                handles.axes33.XLimMode = 'auto';
                handles.axes33.XTickLabelMode = 'auto';
                handles.axes33.XTickLabelRotation=0;
                handles.axes33.YTickMode = 'auto';
                handles.axes33.YLimMode = 'auto';
                handles.axes33.YTickLabelMode = 'auto';
                handles.axes33.XAxis.FontSize = 10;
                handles.axes33.YAxis.FontSize = 10;
                handles.axes33.XAxis.Label.FontSize = 10;
                handles.axes33.YAxis.Label.FontSize = 10;
             end

             handles.axes33.XLim = [0.5 sum(Ix)+0.5]; 
             handles.axes33.YLim = [0.5 sum(Ix)+0.5];
             handles.axes33.XAxis.Label.String = 'Features';
         catch ERR
             errordlg(sprintf('Matrix cannot be displayed!\n%s', ERR.message));
         end
         
    case 'Model P value histogram'

         set(handles.pn3DView,'Visible','off'); set(handles.axes33,'Visible','on'); cla; hold on;
         ah=histogram(handles.axes33,y,'Normalization','probability','EdgeColor','none','FaceColor',rgb('SkyBlue'));%,'FaceAlpha',0.5); 
         maxah= nm_nanmax(ah.Values); ylim([0 maxah]); 
         [f, xi] = ksdensity(y,'function','pdf'); 
         plot(handles.axes33,xi, scaledata(f,[],0, maxah),'Color',rgb('DeepSkyBlue'), 'LineWidth',2);
         handles.axes33.YTick = 0:maxah/10:maxah;
         yticklabels(handles.axes33,'auto')
         [xl,xlb]=nk_GetScaleYAxisLabel(handles.NM.analysis{handles.curranal}.params.TrainParam.SVM);
         xlim(xl); %miny = xl(1); maxy=xl(2);
         xlabel(['Optimization criterion: ' xlb]);
         ylabel('Probability');
         hold on;
         Pval = sum(ve)/size(ve,2);
         xp = [ vp vp ]; yp = [ 0 maxah ];
         if Pval ~= 0
            Pvalstr = sprintf('OOT %s=%1.2f\nOOT Significance: P=%g', xlb, vp, Pval);
         else
            Pvalstr = sprintf('OOT %s=%1.2f\nOOT Significance: P<%g', xlb, vp, 1/perms);
         end
         if exist("vp_cv2","var")
            xp_cv2 = [ vp_cv2 vp_cv2 ];
            patch('Faces', [1 2 3 4], 'Vertices', ...
                [ [ vp_cv2-vp_ci_cv2(1) 0 ]; [ vp_cv2-vp_ci_cv2(1) maxah ] ; [vp_cv2+vp_ci_cv2(1) maxah] ; [vp_cv2+vp_ci_cv2(1) 0 ] ], ...
                'FaceColor', rgb('LightSalmon'), 'EdgeColor', 'none', 'FaceAlpha', 0.3 );
            line( xp_cv2, yp ,'LineWidth',1,'Color',rgb('LightSalmon') );
            Pvalstr_cv2 = sprintf('Mean CV2 %s=%1.2f\nMean CV2 Significance: P=%g', xlb, vp_cv2, ve_cv2);
         end
         line(xp, yp ,'LineWidth',2,'Color','r');
         if exist("vp_cv2","var")
             legend({'Binned null distribution', 'Fitted null distribution', sprintf('95%%-CI Mean %s',xlb), Pvalstr_cv2, Pvalstr }, 'Location', 'northwest');
         else
            legend({'Binned null distribution', 'Fitted null distribution', Pvalstr}, 'Location', 'northwest');
         end

    case NetworkCorrMatStr
         
         set(handles.pn3DView,'Visible','off'); set(handles.axes33,'Visible','on'); cla; hold on
         feats = v.params.features; nF = numel(feats);
         if ~isempty(handles.txtAlterFeats.String)
            altfeatsvar = handles.txtAlterFeats.String;
            feats = evalin('base',altfeatsvar);
         end
       
         sI = 1:nF;
         if filterfl
             FltMetr = handles.selVisMeas2.String{handles.selVisMeas2.Value};
             switch FltMetr
                 case 'CV-ratio of feature weights [Overall Mean]'
                     [ythr,~,~,vlineval] = CVRatio(v,curclass,multiflag);
                     if isempty(filterthr)
                         handles.txtThrVisMeas2.String='-2 2'; 
                         filterthr = vlineval;
                     else
                         filterthr = str2double(filterthr);
                     end
                 case 'CV-ratio of feature weights [Grand Mean]'
                     [ythr,~,~,vlineval] = CVRatio_CV2(v,curclass,multiflag);
                     if isempty(filterthr)
                         handles.txtThrVisMeas2.String='-2 2'; 
                         filterthr = vlineval;
                     else
                         filterthr = str2double(filterthr);
                     end
                 case 'Feature selection probability [Overall Mean]'
                     [ythr,~,~,vlineval] = FeatProb(v,curclass,multiflag);
                     if isempty(filterthr)
                         handles.txtThrVisMeas2.String='0.5'; 
                         filterthr = vlineval;
                     else
                         filterthr = str2double(filterthr);
                     end
                 case 'Probability of feature reliability (95%-CI) [Grand Mean]'
                     [ythr,~,~,vlineval] = Prob_CV2(v,curclass,multiflag);
                     if isempty(filterthr)
                         handles.txtThrVisMeas2.String='-0.5 0.5'; 
                         filterthr = vlineval;
                     else
                         filterthr = str2double(filterthr);
                     end
                 case 'Sign-based consistency -log10(P value)'
                     [ythr,~,~,vlineval] = SignBased_CV2_p_uncorr(v,curclass,multiflag);
                     if isempty(filterthr)
                         handles.txtThrVisMeas2.String=num2str(-log10(0.05),'%1.2f'); 
                         filterthr = vlineval;
                     else
                         filterthr = str2double(filterthr);
                     end
                 case 'Sign-based consistency -log10(P value, FDR)'
                     [ythr,~,~,vlineval] = SignBased_CV2_p_fdr(v,curclass,multiflag);
                     if isempty(filterthr)
                         handles.txtThrVisMeas2.String=num2str(-log10(0.05),'%1.2f'); 
                         filterthr = vlineval;
                     else
                         filterthr = str2double(filterthr);
                     end
                 case 'Analytical -log10(P Value) for Linear SVM [Grand Mean]'
                     [ythr,~,~,vlineval] = Analytical_p(v,curclass,multiflag);
                     if isempty(filterthr)
                         handles.txtThrVisMeas2.String=num2str(-log10(0.05),'%1.2f'); 
                         filterthr = vlineval;
                     else
                         filterthr = str2double(filterthr);
                     end
                 case 'Analytical -log10(P Value, FDR) for Linear SVM [Grand Mean]'
                     [ythr,~,~,vlineval] = Analyitcal_p_fdr(v,curclass,multiflag);
                     if isempty(filterthr)
                         handles.txtThrVisMeas2.String=num2str(-log10(0.05),'%1.2f'); 
                         filterthr = vlineval;
                     else
                         filterthr = str2double(filterthr);
                     end
             end
             if sortfl, [~,sI] = sort(abs(ythr),'descend'); end
             y = y(sI,sI);
             feats = feats(sI);
             y = y(featind,featind);
             feats = feats(featind);
             %y.diag
             if numel(filterthr)>1
                Ix2 = (ythr >= min(filterthr) | ythr <= max(filterthr))';
             else 
                Ix2 = (abs(ythr)>= filterthr)'; 
             end
             Ix2 = Ix2(sI); Ix2 = Ix2(featind);
             Ix = any(y) & Ix2;
         else
             if sortfl, [~,sI] = sort(mean(y),'descend'); end
             y = y(sI,sI);
             feats = feats(sI);
             y = y(featind,featind);
             feats = feats(featind);
             Ix = any(y);
         end
         y = y(Ix,Ix);
         feats = feats(Ix);
         nF = numel(feats);
         if nF<50, FS = 10; elseif nF<100, FS = 7; elseif nF < 150, FS = 5; else, FS = 0; end
         yy = y; 
                        
         try
             cla reset;
             switch meas{measind}
                case NetworkCorrMatStr{1}
                 
                    gy = graph(yy);
                    gy = rmedge(gy, 1:numnodes(gy), 1:numnodes(gy));
                    pgy = plot(gy);
                    pgy.NodeLabel = feats;
                    pgy.EdgeCData = gy.Edges.Weight;
                    pgy.EdgeCData(pgy.EdgeCData > 0) = 0.5;
                    pgy.EdgeCData(pgy.EdgeCData < 0) = -0.5;
                    pgy.EdgeColor = 'flat';
                    colormap gray;
                    gy.Edges.LWidths = 7*abs(gy.Edges.Weight);
                    pgy.LineWidth = gy.Edges.LWidths;
                    pgy.NodeColor = 'k';
                 case NetworkCorrMatStr{2}
                    gy = graph(yy);
                    gy = rmedge(gy, 1:numnodes(gy), 1:numnodes(gy));
                    %nonsigIdx = find(gy.Edges.Weight >= 0.05);
                    %gy = rmedge(gy, nonsigIdx);
                    pgy = plot(gy);
                    pgy.NodeLabel = feats;
                    pgy.EdgeCData = gy.Edges.Weight;
                    pgy.EdgeCData(pgy.EdgeCData >= 0.05) = 0.5;
                    pgy.EdgeCData(pgy.EdgeCData < 0.05) = -0.5;
                    pgy.EdgeAlpha = 0.1;
                    %pgy.EdgeAlpha(pgy.EdgeCData < 0.05) = 0.5;
                    pgy.EdgeColor = 'flat';
                    colormap gray;
                    gy.Edges.LWidths = 0.1*(max(gy.Edges.Weight)+1 -(abs(gy.Edges.Weight)));
                    pgy.LineWidth = gy.Edges.LWidths;
                    pgy.NodeColor = 'k';
                 case NetworkCorrMatStr{3}
                    yy = triu(yy)+triu(yy,1)';
                    gy = graph(upper(yy));
                    gy = rmedge(gy, 1:numnodes(gy), 1:numnodes(gy));
                    %nonsigIdx = find(gy.Edges.Weight >= 0.05);
                    %gy = rmedge(gy, nonsigIdx);
                    pgy = plot(gy);
                    pgy.NodeLabel = feats;
                    pgy.EdgeCData = gy.Edges.Weight;
                    pgy.EdgeCData(pgy.EdgeCData >= 0.05) = 0.5;
                    pgy.EdgeCData(pgy.EdgeCData < 0.05) = -0.5;
                    %pgy.EdgeAlpha = 0.1;
                    %pgy.EdgeAlpha(pgy.EdgeCData < 0.05) = 0.5;
                    pgy.EdgeColor = 'flat';
                    colormap gray;
                    gy.Edges.LWidths = 0.1*(max(gy.Edges.Weight)+1 -(abs(gy.Edges.Weight)));
                    pgy.LineWidth = gy.Edges.LWidths;
                    pgy.NodeColor = 'k';
             end
             %cbar.Label.String=clabel; cbar.Label.FontWeight='bold';


         catch ERR
             errordlg(sprintf('Network plot cannot be displayed!\n%s', ERR.message));
         end
         %handles.axes33.Colormap =  colormap(redgreencmap); 
         
         %handles.axes33.XLim = [min(pgy.XData)-1 max(pgy.XData)+1]; 
         %handles.axes33.YLim = [min(pgy.YData)-1 max(pgy.YData)+1];
        
    otherwise 

        if sortfl 
            [~,ind] = sort(abs(y),'descend');
            y = y(ind); y = y(featind); 
            if exist('se','var'), se = se(ind); se=se(featind); end
        else
            ind = (1:v.params.nfeats)';
            y = y(featind);
            if exist('se','var'), se=se(featind); end
        end

        switch v.params.visflag
            case {0, 3, 4, 5}
                set(handles.pn3DView,'Visible','off'); set(handles.axes33,'Visible','on'); 
                cla
                colorbar('Visible','off')
                handles.axes33.XLimMode = 'auto';
                handles.axes33.XTickLabelMode = 'auto';
                handles.axes33.XTickLabelRotation=0;
                hold on
                if ~isempty(vlineval)
                    y2 = y; 
                    if numel(vlineval)>1
                        y2(y>vlineval(1) & y<vlineval(2))=NaN; 
                    else
                        y2(y<vlineval)=NaN; 
                    end
                    barh(x, y,'FaceColor','b','BarWidth',0.9,'FaceColor',rgb('LightGrey'),'EdgeColor','none','LineWidth',1)
                    barh(x, y2,'FaceColor','b','BarWidth',0.9,'FaceColor',rgb('SkyBlue'),'EdgeColor','none','LineWidth',1)
                else
                    barh(x, y,'FaceColor','b','BarWidth',0.9,'FaceColor',rgb('SkyBlue'),'EdgeColor','none','LineWidth',1)
                end
                if exist('se','var')
                    h = herrorbar(y, x , se, se,'ko');
                    set(h,'MarkerSize',0.001);
                    ylimoff = nm_nanmax(se);
                else 
                    ylimoff = 0;
                end

                if ~isempty(vlineval) 
                    for o=1:numel(vlineval)
                        xline(vlineval(o),'LineWidth',1.5, 'Color', 'red');
                    end
                end

                xlabel(meas{measind},'FontWeight','bold');
                ylabel('Features','FontWeight','bold');
                set(gca,'YTick',x);
                if ~isempty(handles.txtAlterFeats.String)
                    altfeatsvar = handles.txtAlterFeats.String;
                    tfeats = evalin('base',altfeatsvar);
                    if isempty(tfeats)
                        warndlg(sprintf('Could not find the alternative feature vector ''%s'' in the MATLAB workspace!', altfeatsvar))
                        feats = v.params.features;
                    elseif numel(tfeats)~= numel(v.params.features)
                        warndlg(sprintf('The alternative feature vector doesn not have same feature number (n=%g) as the original vector (n=%g)',numel(tfeats),numel(v.params.features))) 
                        feats = v.params.features;
                    elseif ~iscellstr(tfeats) && ~isstring(tfeats)
                        warndlg('The alternative feature vector must be a cell array of strings!'); 
                        feats = v.params.features;
                    else
                        feats = tfeats;
                    end
                else
                    feats = v.params.features;
                end
                
                feats = feats(ind); 
                feats = feats(featind);
                if ~isempty(vlineval)
                    Ix = isnan(y2);
                    featsI = strcat('\color{gray}',feats(Ix));
                    feats(Ix) = featsI;
                end
                feats = regexprep(feats,'_', '\\_');
                handles.axes33.YTickLabel = feats;
                handles.axes33.TickLabelInterpreter = 'tex';
                handles.axes33.YLim = [ x(1)-0.5 x(end)+0.5 ];
                handles.axes33.XLim = [miny-ylimoff maxy+ylimoff ];
                handles.axes33.XTickMode='auto';
                if numel(ind) ~= numel(handles.visdata_table(handles.curmodal,handles.curlabel).tbl.ind)
                    act = 'create';
                else
                    act = 'reorder';
                end
                handles.visdata_table(handles.curlabel, handles.curmodal) = create_visdata_tables(v, handles.visdata_table(handles.curmodal,handles.curlabel), ind, act);

            case 1
                st.fig = handles.pn3DView; 
                st.NMaxes = [ handles.axes26 handles.axes27 handles.axes28];
                set(handles.pn3DView,'Visible','on'); set(handles.axes33,'Visible','off');
                nk_WriteVol(y,'temp',2,v.params.brainmask,[], ...
                    handles.NM.datadescriptor{v.params.varind}.threshval, ...
                    char(handles.NM.datadescriptor{v.params.varind}.threshop));
                if ~isfield(handles,'orthviews')
                    handles.orthviews = nk_orthviews('Image','temp.nii'); 
                else
                    nk_orthviews('Redraw')
                end
                colormap(st.NMaxes(3), jet);
                colormap(st.NMaxes(2), jet);
                colormap(st.NMaxes(1), jet);
                pos = st.NMaxes(3).Position;
                cl = colorbar(st.NMaxes(3));
                cl.TickLabels=linspace(min(y),max(y),6);
                cl.Label.String = meas{measind};
                cl.Label.FontWeight = 'bold';
                cl.Label.FontSize = 11;
                st.NMaxes(3).Position = pos;
                handles.axes33.XLabel.String='';
                handles.axes33.YLabel.String='';
        end


        legend('off');

end

guidata(handles.figure1, handles);

function [y, miny, maxy, vlineval] = CVRatio(v, curclass, multiflag)

if iscell(v.CVRatio)
    if multiflag
        y = nm_nanmean(nk_cellcat(v.CVRatio,[],2),2);
    else
        y = v.CVRatio{curclass}; 
    end
else
    y = v.CVRatio; 
end
vlineval = [-2 2];
miny = nanmin(y); maxy = nanmax(y);

function [y, miny, maxy, vlineval] = CVRatio_CV2(v, curclass, multiflag)

if iscell(v.CVRatio_CV2)
    if multiflag
        y = nm_nanmean(nk_cellcat(v.CVRatio_CV2,[],2),2);
    else
        y = v.CVRatio_CV2{curclass};  
    end
else
    y = v.CVRatio_CV2;
end
vlineval = [-2 2];
miny = nanmin(y); maxy = nanmax(y);

function [y, miny, maxy, vlineval] = FeatProb(v, curclass, multiflag)

if iscell(v.FeatProb)
    if multiflag
        y = nm_nanmean(v.FeatProb{1},2);
    else
        y = v.FeatProb{1}(:,curclass);  
    end   
else
    y = v.FeatProb;
end
vlineval = 0.5;
miny = 0; maxy = nanmax(y);

function [y, miny, maxy, vlineval] = Prob_CV2(v, curclass, multiflag)

if iscell(v.Prob_CV2)
    if multiflag
        y = nm_nanmean(nk_cellcat(v.Prob_CV2,[],2),2);
    else
        y = v.Prob_CV2{curclass};  
    end
else
    y = v.Prob_CV2;
end
vlineval = [-0.5 0.5];
miny = -1; maxy = 1;

function [y, miny, maxy, vlineval] = SignBased_CV2(v, curclass, multiflag)

if iscell(v.SignBased_CV2)
    if multiflag
        y = nm_nanmean(nk_cellcat(v.SignedBased_CV2,[],2),2);
    else
        y = v.SignBased_CV2{curclass};
    end
else
    y = v.SignBased_CV2;
end
vlineval= [];
miny = 0; maxy = 1;

function [y, miny, maxy, vlineval] = SignBased_CV2_z(v, curclass, multiflag)
if iscell(v.SignBased_CV2_z)
    if multiflag
        y = nm_nanmean(nk_cellcat(v.SignBased_CV2_z,[],2),2);
    else
        y = v.SignBased_CV2_z{curclass};
    end
else
    y = v.SignBased_CV2_z;
end
vlineval= [-2 2];
miny = nanmin(y(:)); maxy = nanmax(y(:));

function [y, miny, maxy, vlineval] = SignBased_CV2_p_uncorr(v, curclass, multiflag)

if iscell(v.SignBased_CV2_p_uncorr)
    if multiflag
        y = nm_nanmean(nk_cellcat(v.SignBased_CV2_p_uncorr,[],2),2);
    else
        y = v.SignBased_CV2_p_uncorr{curclass};
    end
else
    y = v.SignBased_CV2_p_uncorr;
end 
miny = 0; maxy = nanmax(y(:));
vlineval = -log10(0.05);

function [y, miny, maxy, vlineval] = SignBased_CV2_p_fdr(v, curclass, multiflag)

if iscell(v.SignBased_CV2_p_fdr)
    if multiflag
        y = nm_nanmean(nk_cellcat(v.SignBased_CV2_p_fdr,[],2),2);
    else
        y = v.SignBased_CV2_p_fdr{curclass};
    end
else
    y = v.SignBased_CV2_p_fdr;
end
miny = 0; maxy = nanmax(y(:));
vlineval = -log10(0.05);

function [y, miny, maxy, vlineval] = Spearman_CV2(v, curclass, multiflag)

if iscell(v.Spearman_CV2)
    if multiflag
        y = nm_nanmean(nk_cellcat(v.Spearman_CV2,[],2),2);
    else
        y = v.Spearman_CV2{curclass};
    end
else
    y = v.Spearman_CV2;
end
miny = nanmin(y(:)); maxy = nanmax(y(:));
vlineval= [];

function [y, miny, maxy, vlineval] = Pearson_CV2(v, curclass, multiflag)

if iscell(v.Pearson_CV2)
    if multiflag
        y = nm_nanmean(nk_cellcat(v.Pearson_CV2,[],2),2);
    else
        y = v.Pearson_CV2{curclass};
    end
else
    y = v.Pearson_CV2;
end
miny = nanmin(y(:)); maxy = nanmax(y(:));
vlineval= [];

function [y, se, miny, maxy, vlineval] = Spearman_CV2_p_uncorr(v, curclass, multiflag)

if iscell(v.Spearman_CV2_p_uncorr)
    if multiflag
        y = nm_nanmean(nk_cellcat(v.Spearman_CV2_p_uncorr,[],2),2);
        if isfield (v,'Spearman_CV2_p_uncorr_STD')
            se = nm_nanmean(sqrt(nk_cellcat(v.Spearman_CV2_p_uncorr_STD,[],2)),2);
        end
    else
        y = v.Spearman_CV2_p_uncorr{curclass}; 
        if isfield (v,'Spearman_CV2_p_uncorr_STD')
            se = sqrt(v.Spearman_CV2_p_uncorr_STD{curclass});
        end
    end
else
    y = v.Spearman_CV2_p_uncorr;
    if isfield (v,'Spearman_CV2_p_uncorr_STD')
        se = sqrt(v.Spearman_CV2_p_uncorr_STD);
    end
end
miny = 0; maxy = nanmax(y(:));
vlineval = -log10(0.05);

function [y, se, miny, maxy, vlineval] = Pearson_CV2_p_uncorr(v, curclass, multiflag)

if iscell(v.Pearson_CV2_p_uncorr)
    if multiflag
        y = nm_nanmean(nk_cellcat(v.Pearson_CV2_p_uncorr,[],2),2);
        if isfield (v,'Pearson_CV2_p_uncorr_STD')
            se = nm_nanmean(sqrt(nk_cellcat(v.Pearson_CV2_p_uncorr_STD,[],2)),2);
        end
    else
        y = v.Pearson_CV2_p_uncorr{curclass};  
        if isfield (v,'Pearson_CV2_p_uncorr_STD')
            se = sqrt(v.Pearson_CV2_p_uncorr_STD{curclass});
        end
    end
else
    y = v.Pearson_CV2_p_uncorr;
    if isfield (v,'Pearson_CV2_p_uncorr_STD')
        se = sqrt(v.Pearson_CV2_p_uncorr_STD);
    end
end
miny = 0; maxy = nanmax(y(:));
vlineval = -log10(0.05);

function [y, miny, maxy, vlineval] = Spearman_CV2_p_fdr(v, curclass, multiflag)
if iscell(v.Spearman_CV2_p_fdr)
    if multiflag
        y = nm_nanmean(nk_cellcat(v.Spearman_CV2_p_fdr,[],2),2);
    else
        y = v.Spearman_CV2_p_fdr{curclass}; 
    end
else
    y = v.Spearman_CV2_p_fdr;
end
miny = 0; maxy = nanmax(y(:));
vlineval = -log10(0.05);

function [y, miny, maxy, vlineval] = Pearson_CV2_p_fdr(v, curclass, multiflag)
if iscell(v.Pearson_CV2_p_fdr)
    if multiflag
        y = nm_nanmean(nk_cellcat(v.Pearson_CV2_p_fdr,[],2),2);
    else
        y = v.Pearson_CV2_p_fdr{curclass};  
    end
else
    y = v.Pearson_CV2_p_fdr;
end
miny = 0; maxy = nanmax(y(:));
vlineval = -log10(0.05);

function [y, miny, maxy, vlineval] = PermZ_CV2(v, curclass, multiflag)
if iscell(v.PermZ_CV2)
    if multiflag
        y = nm_nanmean(nk_cellcat(v.PermZ_CV2,[],2),2);
    else
        y = v.PermZ_CV2{curclass};  
    end
else
    y = v.PermZ_CV2;
end
miny = nanmin(y(:)); maxy = nanmax(y(:));
vlineval = [];

function [y, miny, maxy, vlineval] = PermProb_CV2(v, curclass, multiflag)
if iscell(v.PermProb_CV2)
    if multiflag
        y = nm_nanmean(nk_cellcat(v.PermProb_CV2,[],2),2);
    else
        y = v.PermProb_CV2{curclass};  
    end
else
    y = v.PermProb_CV2;
end
miny = 0; maxy = nanmax(y(:));
vlineval = -log10(0.05);

function [y, miny, maxy, vlineval] = PermProb_CV2_FDR(v, curclass, multiflag)
if iscell(v.PermProb_CV2_FDR)
    if multiflag
        y = nm_nanmean(nk_cellcat(v.PermProb_CV2_FDR,[],2),2);
    else
        y = v.PermProb_CV2_FDR{curclass};  
    end
else
    y = v.PermProb_CV2_FDR;
end
miny = 0; maxy = nanmax(y(:));
vlineval = -log10(0.05);

function [y, miny, maxy, vlineval] = Analytical_p(v, curclass, multiflag)
if iscell(v.Analytical_p)
    if multiflag
        y = nm_nanmean(nk_cellcat(v.Analytical_p,[],2),2);
    else
        y = v.Analytical_p{curclass};  
    end
else
    y = v.Analytical_p;
end
miny = 0; maxy = nanmax(y(:));
vlineval = -log10(0.05);

function [y, miny, maxy, vlineval] = Analyitcal_p_fdr(v, curclass, multiflag)
if iscell(v.Analyitcal_p_fdr)
    if multiflag
        y = nm_nanmean(nk_cellcat(v.Analyitcal_p_fdr,[],2),2);
    else
        y = v.Analyitcal_p_fdr{curclass};  
    end
else
    y = v.Analyitcal_p_fdr;
end
miny = 0; maxy = nanmax(y(:));
vlineval = -log10(0.05);

function [y, se, miny, maxy, vlineval] = MEAN(v, curclass, multiflag)
if iscell(v.MEAN)
    if multiflag
        y = nm_nanmean(nk_cellcat(v.MEAN,[],2),2);
        se = nm_nanmean(nk_cellcat(v.SE,[],2),2);
    else
        y = v.MEAN{curclass}; se = v.SE{curclass}; 
    end
else
    y = v.MEAN; se = v.SE; 
end
miny = nanmin(y(:)); maxy = nanmax(y(:));
vlineval = [];

function [y, se, miny, maxy, vlineval] = MEAN_CV2(v, curclass, multiflag)
if iscell(v.MEAN_CV2)
    if multiflag
        y = nm_nanmean(nk_cellcat(v.MEAN_CV2,[],2),2);
        se = nm_nanmean(nk_cellcat(v.SE_CV2,[],2),2);
    else
        y = v.MEAN_CV2{curclass}; se = v.SE_CV2{curclass}; 
    end
else
    y = v.MEAN_CV2; se = v.SE_CV2; 
end
miny = nanmin(y(:)); maxy = nanmax(y(:));
vlineval = [];