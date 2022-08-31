function imagingMLI(casenum, MLIcont, APPstMLI, brainmask, badcoords, datadescriptor)
global stMLI

y = MLIcont.Y_mapped(casenum,:)';
stMLI = APPstMLI;
nk_WriteVol(y,'tempMLI', 2, brainmask, badcoords, datadescriptor.threshval, char(datadescriptor.threshop),[],[],false);
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