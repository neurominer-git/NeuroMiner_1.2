function NM = nk_AdjBrainMask(NM)

i = nk_input('Define modality index',0,'i');
if i>numel(NM.Y), error('Index exceeds no. of data modalities'); end
[file_brainmask,path_brainmask] = uigetfile({'*.nii';'*.img';'*.*'},...
                          'Brainmask Image Selector');
[~,n,e] = fileparts(NM.datadescriptor{i}.input_settings.brainmask);
old_file = [n e];
if ~strcmp(file_brainmask, old_file)
    warning('name of original brainmask image and new weighting image is not identical.')
end
pth = fullfile(path_brainmask,file_brainmask);
NM.brainmask{i} = pth;
NM.datadescriptor{i}.input_settings.brainmask = pth;
NM.datadescriptor{i}.input_settings.Vm = spm_vol(pth);
if isfield(NM.datadescriptor{i}.input_settings,'Pw') && ~isempty(NM.datadescriptor{i}.input_settings.Pw)
    [p,n,e] = fileparts(NM.datadescriptor{i}.input_settings.Pw);
    old_file = [n e];
    [file_wmask,path_wmask] = uigetfile({'*.nii';'*.img';'*.*'},...
                          'Weighting Image Selector');
    if ~strcmp(file_wmask, old_file)
        warning('name of original weighting image and new weighting image is not identical.')
    end
    pth = fullfile(path_wmask,file_wmask);
    NM.datadescriptor{i}.input_settings.Pw = pth;
    NM.datadescriptor{i}.input_settings.Vw = spm_vol(pth);
end

