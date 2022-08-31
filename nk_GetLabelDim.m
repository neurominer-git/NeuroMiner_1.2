function nl = nk_GetLabelDim(MULTILABEL)

if isfield(MULTILABEL,'sel')
    nl = numel(MULTILABEL.sel);
else
    nl = MULTILABEL.dim;
end