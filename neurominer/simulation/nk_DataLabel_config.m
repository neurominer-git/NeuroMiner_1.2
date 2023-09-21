function altlabels = nk_DataLabel_config(n_subjects_all, altlabelnames)
% this function asks the user to define alternative labels for the
% validation data as used in the locked analyses

unique_altlabels = unique(altlabelnames);
c = cell(length(unique_altlabels),1);
altlabels = cell2struct(c, unique_altlabels);


for i=1:length(unique_altlabels)
    promptstr = sprintf('Enter %s labels for validation sample', unique_altlabels{i});
    altlabels.(unique_altlabels{i}) = nk_input(promptstr,0,'r',[],[n_subjects_all 1]);
end
