function res = JuSpace_noGUI_2D(list1,atlas,options,PET_list,image_for_size,dir_tool)


% addpath(genpath('/opt/matlab/JuSpace_v1.3/'));

list2 = [];

if ~ischar(image_for_size)
   image_for_size = char(image_for_size); 
end

% list1: cell with filenames
% list2 (options

% JuSpace directory
% dir_tool = fileparts(which('JuSpace'));

% atlas file
% atlas = fullfile(dir_tool,'atlas',atlas);

% PET directory and selected files
dir_PET = fullfile(dir_tool,'PETatlas/');
aa = dir(dir_PET);
names_PET = {aa(3:end).name}';
files_PET_all = strcat(repmat(dir_PET, size(aa(3:end),1),1), names_PET);
filesPET = files_PET_all(PET_list,1);


for i = 1:size(names_PET,1)
   tt = regexp(names_PET{i},'_','Split');
   Rec_list{i}=tt{1};
end

PET_names = Rec_list(PET_list);

opt_comp = options(1);
opt_ana = options(2);
opt_perm = options(3);
opt_auto = options(4);
opt_perm_spatial = options(5);

file_part = '';
switch options(1)
    case 1
        file_part = [file_part '_ESb'];
    case 2
        file_part = [file_part '_ESw'];
    case 3
        file_part = [file_part '_mList1'];
    case 4
        file_part = [file_part '_List1Each'];
    case 5
        file_part = [file_part '_indZ'];
    case 6
        file_part = [file_part '_pwDiff'];
    case 7
        file_part = [file_part '_looList1'];
    case 8
        file_part = [file_part '_List1EachGroupTest'];
end

switch options(2)
    case 1
        file_part = [file_part '_Spearman'];
    case 2
        file_part = [file_part '_Pearson'];
    case 3
        file_part = [file_part '_multReg'];
end

if options(3) == 1 
    file_part = [file_part '_withExactP'];
end

if options(5) == 1
    file_part = [file_part '_withExactSpatialP'];
end

% global files_PET
% aa = get(handles.PETlist,'Value');
% str_PET = get(handles.PETlist,'String');
% str_PET = str_PET(aa);
% for i = 1:length(str_PET);
%    tt = regexp(str_PET{i},'_','Split');
%    Rec_list{i}=tt{1};
% end
% filesPET = files_PET(aa);
% atlas = char(get(handles.atlaslist,'String'));

time_now = datestr(datetime('now'),'ddmmmyyyy_HHMMSS');

% file_save = fullfile(dir_save,['Results_',name_save file_part '_' time_now '.mat']);

if options(1)<4
%     image_save = fullfile(dir_save,[name_save file_part '_' time_now '.nii']);
    [res,p_all,stats,data,D1,D2, data_PET,Resh,T1] = compute_DomainGauges_2D(list1,list2,filesPET,atlas,options,image_for_size,dir_tool,image_save);
else
    [res,p_all,stats,data,D1,D2, data_PET,Resh,T1] = compute_DomainGauges_2D(list1,list2,filesPET,atlas,options,image_for_size,dir_tool);
end

% opt_for_perm = [1,2,5,6];
% opt_for_spat_perm = [3,4,8];
% 
% if options(3)==1 && ismember(options(1),opt_for_perm)% && options(2)~=3
%     disp('Computing exact p-value');
%    [p_exact,dist_r] = compute_exact_pvalue(D1,D2,data_PET,res,Nperm,options,T1);
%    Resh(:,end+1) = [{'p_exact'}; num2cell_my(p_exact')];
% end
% 
% if options(5)==1 && ismember(options(1),opt_for_spat_perm)
%     disp('Computing exact spatial p-value')
%     [p_exact,dist_r] = compute_exact_spatial_pvalue(D1,data_PET,atlas,res,Nperm,options,filesPET, T1);
%     Resh(:,end+1) = [{'p_exact'}; num2cell_my(p_exact')];
% end


end
