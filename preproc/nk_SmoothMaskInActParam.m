function sYw = nk_SmoothMaskInActParam(uPREPROC, PV)

I = arrayfun( @(j) isfield(uPREPROC.ACTPARAM{j},'RANK'), 1:numel( uPREPROC.ACTPARAM ));
sYw=[];
if any(I)
    Ix = find(I);
    for qx = 1:numel(Ix)
        if isfield(uPREPROC.ACTPARAM{Ix(qx)}.RANK,'EXTERN')
            sYw = nk_PerfSpatFilt( uPREPROC.ACTPARAM{Ix(qx)}.RANK.EXTERN, uPREPROC, PV ); 
            %here, we assume that there is only one
            %weighting map to be smoothed alongside the
            %data. This will obviously not work for
            %multiple weighting maps...
            break
        end
    end
end