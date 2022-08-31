function [ inp ] = nk_ApplyLabelTransform( PREPROC, MODEFL, inp )
global MULTILABEL
if MULTILABEL.flag
    if isfield(MULTILABEL,'sel')
        lb=MULTILABEL.sel(inp.curlabel);
    else
        lb=inp.curlabel;
    end
else
    lb = 1;
end

[ inp.label, inp.targscale, inp.minLbCV, inp.maxLbCV, ~, inp.PolyFact ] = nk_LabelTransform(PREPROC, MODEFL, inp.labels(:,lb));

if isfield(inp,'labelOOCV')
    inp.curlabelOOCV = inp.labelOOCV(:,lb);
    if inp.targscale
        IN.minY = inp.minLbCV; IN.maxY = inp.maxLbCV; 
        inp.curlabelOOCV = nk_PerfScaleObj(inp.labelOOCV(:, lb), IN); 
    end
    if ~isempty(inp.PolyFact), inp.curlabelOOCV = inp.labelOOCV(:,inp.curlabel) .^ (1/inp.PolyFact); end 
end



