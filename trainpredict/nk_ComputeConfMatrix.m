function confmatrix = nk_ComputeConfMatrix(label, pred, ngroups)

ind = ~isnan(pred) & ~isnan(label);
confmatrix = zeros(ngroups);
lx = size(label,1);
% Compute confusion matrix
u=unique(label);
for i=1:lx
    if ~ind(i), continue, end
    iu = u==label(i);
    ip = u==pred(i);
    confmatrix(iu,ip) = confmatrix(iu,ip) + 1;
end
   
