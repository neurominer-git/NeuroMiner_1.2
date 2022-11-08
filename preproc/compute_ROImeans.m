function Ymean = compute_ROImeans(Yimg, brainmask, atlas)%, PREPROC, prevP)

S.Vm                         = spm_vol(brainmask);
[S.dims, S.indvol, ~, S.vox] = nk_ReadMaskIndVol(S.Vm, []);

image_for_size = char(brainmask);
atlas = [];
if ~exist('atlas','var') || isempty(atlas)
    maskVec = nk_Vol2Vec(brainmask);
else 
    ROI_matrix = resize_img_useTemp_imcalc_NM(char(atlas),image_for_size); %from JuSpace Toolbox
    maskVec = reshape(ROI_matrix,size(ROI_matrix,1)*size(ROI_matrix,2)*size(ROI_matrix,3),1)';
end 
thres = 0; % should be possible to configure in menu
maskVec = int64(maskVec(maskVec > thres));

% project Yimg back to input space 
%bpYimg = project2InputSpace(Yimg, PREPROC, prevP);

% create mean GMV table
Ymean = zeros(size(Yimg,1),int64(max(maskVec)));
for i = 1:max(maskVec)
    indVec = maskVec == i;
    Ymean(:,i) = mean(Yimg(:,indVec),2);
end


end