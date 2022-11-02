function Ypet = apply_JuSpace(Yimg, brainmask, atlas, cortype, autocorcorrect, petlist, dir_save)
global JSMEM 

S.Vm                         = spm_vol(brainmask);
[S.dims, S.indvol, ~, S.vox] = nk_ReadMaskIndVol(S.Vm, []);
image_for_size = string(brainmask);

if ~isfield(JSMEM,'ROI_matrix') || isempty(JSMEM.ROI_matrix)
    ROI_matrix = resize_img_useTemp_imcalc_NM(char(atlas),char(image_for_size),dir_save); %from JuSpace Toolbox
    JSMEM.ROI_matrix = ROI_matrix;
else
    ROI_matrix = JSMEM.ROI_matrix;
end 
atlasVec = reshape(ROI_matrix,size(ROI_matrix,1)*size(ROI_matrix,2)*size(ROI_matrix,3),1)';

V = zeros(S.dims);
C = zeros(size(Yimg,1),numel(V));
for i = 1:size(Yimg,1)
    V = zeros(S.dims);
    V(S.indvol) = Yimg(i,:); % transfer vector into 3D space
    C(i,:) = reshape(V,1,numel(V));
end

a = unique(atlasVec(:));
a = a(a~=0);
atlas_vals = a(~isnan(a));

% create mean GMV table
Ymean = zeros(size(Yimg,1),numel(atlas_vals));
clear Yimg;
clear V;
clear S;
clear ROI_matrix;
for i = 1:numel(atlas_vals)
    indVec = round(atlasVec) == round(atlas_vals(i));
    Ymean(:,i) = mean(C(:,indVec),2);
end

% save_dir = NM.
% V = zeros(S.dims);
% C = zeros(size(Yimg,1),numel(V));
% V = nk_WriteVol(Yimg, 
% for i = 1:size(Yimg,1)
%     % create temp directory 
%     if ~exist('tempDir','dir')
%         mkdir('tempDir');
%     end
%     
%     % make temp filenames 
%     
%     nk_WriteVol(Yimg, 'tempDir
%     V = zeros(S.dims);
%     V(S.indvol) = Yimg(i,:); % transfer vector into 3D space
%     C(i,:) = reshape(V,1,numel(V));
% end
 
options = [4; cortype; 0; autocorcorrect; 0]; 

petvec = zeros([1,numel(petlist)]);
for i = 1:numel(petlist)
    petvec(i) = petlist{i}.listidx; 
end

Ypet = JuSpace_noGUI_2D(Ymean,atlas,options,petvec,image_for_size,dir_save);

    
end
