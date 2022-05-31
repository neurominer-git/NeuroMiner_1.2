function STATUS = nk_CheckFieldStatus(parent, children, grandchildren, grandchildrencrit, STATUS)

strdef = '...'; strundef = '???'; strcrit = '!!!';

n = numel(grandchildren);
missdef = cellstr(repmat(strundef,n,1));
okdef   = cellstr(repmat(strdef,n,1));
if exist('grandchildren','var') && ~isempty(grandchildren)
    if exist('grandchildrencrit','var') && ~isempty(grandchildrencrit)
        missdef{grandchildrencrit} = strcrit;
    end
else
    grandchildren = [];
end

if isempty(parent)
    childrenfields = [];
    nC = 0;
else
    childrenfields = fieldnames(parent);
    nC = numel(childrenfields);
end

newfl = 1;
for i=1:numel(children)
    STATUS.(children{i}) = strundef;
    for j=1:nC
        if strcmp(childrenfields{j},children{i}) 
            newfl = 0; STATUS.(children{i}) = strdef; break
        end
    end
end

if ~isempty(grandchildren) && nC >0
    grandchildrenfields = fieldnames(parent.(children{1}));
    for i=1:numel(grandchildren)
        if newfl
            STATUS.(grandchildren{i}) = missdef{i};
        else
            fl = 1;
            for j=1:numel(grandchildrenfields)
                if strcmp(grandchildrenfields{j},grandchildren{i}); fl = 0; break; end
            end
            if fl
                STATUS.(grandchildren{i}) = missdef{i};
            else
                STATUS.(grandchildren{i}) = okdef{i};
            end
        end
    end
end

