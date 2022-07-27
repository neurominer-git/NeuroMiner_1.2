function imagingMLI(casenum, MLIcont, APPstMLI)
%global stMLI
global stMLI
stMLI = APPstMLI;
y = MLIcont.Y_mapped(casenum,:)';
%y_ciu = MLIcont.Y_mapped_ciu(casenum,:);
%y_cil = MLIcont.Y_mapped_cil(casenum,:);

brainmask = evalin('base', 'NM.brainmask{1}'); % curModal 

nk_WriteVol(y,'temp', 2, brainmask, [], 0, 'gt');
mli_orthviews('Image','temp.nii'); 


end