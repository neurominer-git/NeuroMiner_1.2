function rT = ranktransform(tr, ts)

if ~exist('ts','var') || isempty(ts)
    ts = tr;
end
rT = zeros(size(ts));
for i=1:size(rT,1)
    rT(i,:) = nk_ComputePercentiles(tr,ts(i,:),'inverse');
end