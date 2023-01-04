function NM = update_NMstruct(NM)
% function that updates an older NM structure to be compatible with the new
% release

% add label-field to analysis if necessary

% create default label struct
label.label = NM.labels; 
label.modeflag = NM.modeflag; 
label.altlabelflag = 0; 

for i=1:length(NM.analysis)
    if ~isfield(NM.analysis{1,i}.params, 'label')
        NM.analysis{1,i}.params.label = label;
    end
end

