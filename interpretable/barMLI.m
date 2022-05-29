function [f, ax] = barMLI(casenum, MLIcont, feats)

f=figure; hold on; 
nF = numel(feats);
y = MLIcont.Y_mapped(casenum,:);
y_ciu = MLIcont.Y_mapped_ciu(casenum,:);
y_cil = MLIcont.Y_mapped_cil(casenum,:);
[~,idx] = sort(y);
idx2 = abs(y_ciu) > abs(y);
if any(idx2)
    y2 = y;
    y2(~idx2) = nan;
    y(idx2) = nan;
end

bar(y(idx),'FaceColor',rgb('DodgerBlue'),'EdgeColor',rgb("SlateGrey")); 
errorbar(1:nF, y(idx), y_ciu(idx), y_cil(idx), 'LineStyle','none', 'Color', 'k'); 
if any(idx2)
    bar(y2(idx), 'FaceColor',rgb('LightGrey'),'EdgeColor',rgb("Grey")); 
    errorbar(1:nF, y2(idx), y_ciu(idx), y_cil(idx), 'LineStyle','none', 'Color', rgb("Grey")); 
end
ax=gca; 
ax.XTick=1:numel(feats); 
if any(idx2)
    featsI = strcat('\color{gray}',feats(idx2));
    feats(idx2) = featsI;
end
feats = regexprep(feats,'_', '\\_');
ax.XTickLabel=feats(idx);
ax.Title.String = sprintf('Predictive profile of subject #%g', casenum);
ax.Title.FontWeight = 'bold';
ax.YAxis.Label.String = 'Percentage change of prediction [ % maximum range ]';
ax.XAxis.Label.String = 'Features'; 
ax.Box="on";
f.Position = [100 100 750 750];