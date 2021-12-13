function LabelInteraction = nk_Check4LabelInteraction(P)
LabelInteraction = false; 
for a = 1:numel(P)
    if P{a}.LabelInteraction, LabelInteraction = true; end
end