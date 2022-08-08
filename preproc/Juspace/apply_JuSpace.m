function Ypet = apply_JuSpace(Yimg, brainmask, atlas, cortype, autocorcorrect, petlist)

S.brainmask                  = brainmask;
S.Vm                         = spm_vol(S.brainmask);

[S.dims, S.indvol, ~, S.vox] = nk_ReadMaskIndVol(S.Vm, []);
% for now, you cannot prune any columns before! (add later --> similar to
% smoothing function

ttY                      = Yimg;

%nanMask                      = isnan(ttY);preprpr
%indnan                       = any(nanMask,2);

V = zeros(S.dims);
C = zeros(size(ttY,1),numel(V));
for i = 1:size(ttY,1)
    V = zeros(S.dims);
    V(S.indvol) = ttY(i,:); % transfer vector into 3D space
    C(i,:) = reshape(V,1,numel(V));
end
 
options = [4; cortype; 0; autocorcorrect; 0]; 

petvec = zeros([1,numel(petlist)]);
for i = 1:numel(petlist)
    petvec(i) = petlist{i}.listidx; 
end

image_for_size = string(S.brainmask);
%image_for_size = '/volume/projects/LH_BEST/Data/BEST/MRI/ResCov/ResCov_mwp1368_MRI_sMRI_400096.nii';
Ypet = JuSpace_noGUI_2D(C,atlas,options,petvec,image_for_size);
clear ttY
clear C
    
end