function imagingMLI(casenum, MLIcont, APPstMLI)

global stMLI
stMLI = APPstMLI;
y = MLIcont.Y_mapped(casenum,:)';
%y_ciu = MLIcont.Y_mapped_ciu(casenum,:);
%y_cil = MLIcont.Y_mapped_cil(casenum,:);

brainmask = evalin('base', 'NM.brainmask{1}'); % curModal 

nk_WriteVol(y,'tempMLI', 2, brainmask, [], 0, 'gt');
mli_orthviews('Image','tempMLI.nii'); 
% alphamask = y ~= 0;
% hnd = imagesc(y);
% set(hnd,'AlphaData',alphamask);
colormap(stMLI.NMaxes(3), jet);
colormap(stMLI.NMaxes(2), jet);
colormap(stMLI.NMaxes(1), jet);
pos = stMLI.NMaxes(3).Position;
cl = colorbar(stMLI.NMaxes(3));
%clim([min(y),max(y)]);
cl.TickLabels=round(linspace(min(y),max(y),6),2);
%cl.Label.String = meas{measind};
cl.Label.FontWeight = 'bold';
cl.Label.FontSize = 10;
stMLI.NMaxes(3).Position = pos;
% handles.axes33.XLabel.String='';
% handles.axes33.YLabel.String='';

end