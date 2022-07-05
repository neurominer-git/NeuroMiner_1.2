function [f, ax] = barMLI(casenum, MLIcont, feats, refdata, sortfl, uiaxes, uiaxes2) %newdata

if exist('uiaxes','var')
    ax = uiaxes;
    f=0;
else
    f=figure;
    f.Position = [100 100 750 750];
    ax = gca;
end

hold(ax,'on');
nF = numel(feats);
y = MLIcont.Y_mapped(casenum,:);
y_ciu = MLIcont.Y_mapped_ciu(casenum,:);
y_cil = MLIcont.Y_mapped_cil(casenum,:);
if sortfl
    [~,idx] = sort(y);
else
    idx = 1:size(y,2);
end
idx2 = abs(y_ciu) > abs(y);
if any(idx2)
    y2 = y;
    y2(~idx2) = nan;
    y(idx2) = nan;
end

b1 = bar(ax, y(idx),'FaceColor',rgb('DodgerBlue'),'EdgeColor',rgb("SlateGrey"));
e1 = errorbar(ax,1:nF, y(idx), y_ciu(idx), y_cil(idx), 'LineStyle','none', 'Color', 'k');

legend(ax, b1, 'significant prediction change');
if any(idx2)
    b2 = bar(ax,y2(idx), 'FaceColor',rgb('LightGrey'),'EdgeColor',rgb("Grey"));
    e2= errorbar(ax,1:nF, y2(idx), y_ciu(idx), y_cil(idx), 'LineStyle','none', 'Color', rgb("Grey"));
    legend(ax, [b1,b2], 'significant prediction change', 'non-significant prediction change');
end


%ax=gca;
ax.XTick=1:numel(feats);
if any(idx2)
    featsI = strcat('\color{gray}',feats(idx2));
    feats(idx2) = featsI;
end
feats = regexprep(feats,'_', '\\_');
ax.XTickLabel=feats(idx);
ax.YAxis.Label.String = '95% CI Mean Percentage Change of Prediction';%'Percentage change of prediction [ % maximum range ]';
ax.YAxis.Label.FontWeight = 'bold';
ax.YAxis.Label.FontSize = 10;
ax.XAxis.Label.String = 'Features';
ax.Box="on";

if exist('refdata','var') && ~isempty(refdata)
    if exist('uiaxes2','var')
        ax2 = uiaxes2;
    else
        ax.Position(4) = .50;
        pos = ax.Position;
        pos([2 4]) = [0.75 0.2];
        ax2 = axes('Position', pos);
    end
    if exist('origdata','var')&& ~isempty(origdata) % for CARE project
        refdata = [origdata(casenum,:), refdata];
        centiles = nk_ComputePercentiles(refdata,origdata(casenum,:),'inverse');
    else
        centiles = nk_ComputePercentiles(refdata,refdata(casenum,:),'inverse');
        grouplabel = evalin('base', 'NM.label');
        group1_idx = grouplabel == 1;
        group2_idx = grouplabel == 2; 
        group1 = refdata(group1_idx,:);
        for i = 1:size(group1,1)
            centiles_i = nk_ComputePercentiles(refdata, group1(i,:), 'inverse');
            if i == 1
                centiles_group1 = centiles_i;
            else 
                centiles_group1 = [centiles_group1; centiles_i];
            end
        end
        av_centiles_group1 = mean(centiles_group1);

        group2 = refdata(group2_idx,:);
        for i = 1:size(group2,1)
            centiles_i = nk_ComputePercentiles(refdata, group2(i,:), 'inverse');
            if i == 1
                centiles_group2 = centiles_i;
            else 
                centiles_group2 = [centiles_group2; centiles_i];
            end
        end
        av_centiles_group2 = mean(centiles_group2);
        
    end
    
    bar(ax2, centiles(idx),'FaceColor',rgb('SlateGray'),'EdgeColor',rgb("Black"));
    hold(ax2,'on');
    plot(ax2, idx, av_centiles_group1(idx), 'LineWidth', 3, 'color', rgb('blue'));%,'FaceColor',rgb('SlateGray'),'EdgeColor',rgb("Black"));
    plot(ax2, idx, av_centiles_group2(idx), 'LineWidth', 3, 'color', rgb('orange'));
    hold(ax2, 'off');

    groupnames = evalin('base', 'NM.groupnames');
    legend(ax2, "subject's percentiles" , sprintf('mean percentiles of group %s', char(groupnames(1))), ...
        sprintf('mean percentiles of group %s', char(groupnames(2))));

    ax2.YAxis.Label.String = 'Percentile rank [%]';
    ax2.YAxis.Label.FontWeight = 'bold';
    ax2.YAxis.Label.FontSize = 10;
    ax2.XTick=1:numel(feats);
    ax2.XLim=[-0.2 numel(feats)+1.2];
    ax2.XTickLabel=[];
    if ~exist('uiaxes2','var')
        g = grouplabel(casenum);
        case_groupname = char(groupnames(g));
        ax2.Title.String = sprintf('Predictive profile of subject #%g, Group: %s', casenum, case_groupname);
        ax2.Title.FontWeight = 'bold';
        ax2.Title.FontSize = 14;
    end
else
    if ~exist('uiaxes','var')
        ax.Title.String = sprintf('Predictive profile of subject #%g', casenum);
    end
    ax.Title.FontWeight = 'bold';
    ax.Title.FontSize = 14;
end
