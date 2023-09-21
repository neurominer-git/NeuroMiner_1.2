function [paramstr, paramcell] = print_struct_params(st, delim)

f = fieldnames(st);
nF = numel(f);
paramcell = cell(1,nF);
for i=1:nF
    val = st.(f{i});
    if isstruct(val)
        paramcell{i} = sprintf('%s: struct (%g param)', f{i}, numel(fieldnames(val)));
    elseif ~isfinite(val) || isempty(val)
        paramcell{i} = sprintf('%s: undefined', f{i});
    elseif isnumeric(val)
        if numel(val)>1
            paramcell{i} = sprintf('%s: %s', f{i}, nk_ConcatParamstr(val)); 
        else
            paramcell{i} = sprintf('%s: %g', f{i}, val);
        end
    elseif islogical(val)
        valstr = {'true','false'};
        paramcell{i} = sprintf('%s: %s', f{i}, valstr{val});
    elseif isstring(val) || ischar(val)
        paramcell{i} = sprintf('%s: %s', f{i}, val);
    end
end
paramstr = strjoin(paramcell, delim);