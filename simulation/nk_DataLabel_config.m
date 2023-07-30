function altlabels = nk_DataLabel_config(n_subjects_all, altlabelnames)
% this function asks the user to define alternative labels for the
% validation data as used in the locked analyses

c = cell(length(fields),1);
altlabels = cell2struct(c,altlabelnames);

for i=1:length(altlabelnames)
    promptstr = sprintf('Enter %s labels for validation sample', altlabelnames(i));
    altlabels.(altlabelname(i)) = nk_input(promptstr,0,'r',[],[n_subjects_all 1]);
end
