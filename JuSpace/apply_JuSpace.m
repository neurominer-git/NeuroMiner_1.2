function Ypet = apply_JuSpace(Yimg, brainmask, atlas, cortype, autocorcorrect, petlist)

% add JuSpace path
dir_tool = '/opt/matlab/JuSpace_v1.3/';

% create global variable to store atlas matrix
global JSMEM

% create atlas matrix if not already stored in JSMEM
if ~isfield(JSMEM,'atlas_matrix') || isempty(JSMEM.atlas_matrix)
    atlas_matrix = resize_image_JuSpace(atlas,brainmask);    
    JSMEM.atlas_matrix = atlas_matrix;
else
    atlas_matrix = JSMEM.atlas_matrix;
end

% reshape atlas into 1-D format
atlasVec = reshape(atlas_matrix,size(atlas_matrix,1)*size(atlas_matrix,2)*size(atlas_matrix,3),1)';

% read brainmask for resizing of Yimg
S.Vm                         = spm_vol(brainmask);
[S.dims, S.indvol, ~, S.vox] = nk_ReadMaskIndVol(S.Vm, []);

% resize Yimg to to brainmask
V = zeros(S.dims);
C = zeros(size(Yimg,1),numel(V));
for i = 1:size(Yimg,1)
    V = zeros(S.dims);
    V(S.indvol) = Yimg(i,:); % transfer vector into 3D space
    C(i,:) = reshape(V,1,numel(V));
end

% get number of atlas regions
a = unique(atlasVec(:));
a = a(a~=0);
atlas_vals = a(~isnan(a));

%initialize Ymean with zeros
Ymean = zeros(size(Yimg,1),numel(atlas_vals));

clear Yimg;
clear V;
clear S;
clear atlas_matrix;

for i = 1:numel(atlas_vals)
    indVec = round(atlasVec) == round(atlas_vals(i));
    Ymean(:,i) = mean(C(:,indVec),2);
end
 
options = [4; cortype; 0; autocorcorrect; 0];

petvec = zeros([1,numel(petlist)]);
for i = 1:numel(petlist)
    petvec(i) = petlist{i}.listidx; 
end

Ypet = JuSpace_noGUI_2D(Ymean,atlas,options,petvec,brainmask,dir_tool);
    
end