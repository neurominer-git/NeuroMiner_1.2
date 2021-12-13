function [ inp ] = nk_ApplyLabelTransform( PREPROC, MODEFL, inp )

[ inp.label, inp.targscale, inp.minLbCV, inp.maxLbCV, ~, inp.PolyFact ] = nk_LabelTransform(PREPROC, MODEFL, inp.labels(:,inp.curlabel));
if isfield(inp,'labelOOCV')
    inp.curlabelOOCV = inp.labelOOCV(:,inp.curlabel);
    if inp.targscale, 
        IN.minY = inp.minLbCV; IN.maxY = inp.maxLbCV; 
        inp.curlabelOOCV = nk_PerfScaleObj(inp.labelOOCV(:, inp.curlabel), IN); 
    end
    if ~isempty(inp.PolyFact), inp.curlabelOOCV = inp.labelOOCV(:,inp.curlabel) .^ (1/inp.PolyFact); end 
end



