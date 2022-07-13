function centiles = nk_ComputePercentiles(refdata, values, direction, groupnames)
% =========================================================================
% FORMAT centiles = nk_ComputePercentiles(data, values)
% =========================================================================
% computes the percentile
nd = size(refdata,2);
centiles = zeros(1,nd);
switch direction
    case 'normal'
        for i=1:nd
            centiles(i) = prctile(refdata(:,i), values(i));
        end
    case 'inverse'
        for i=1:nd
            centiles(i) = invprctile(refdata(:,i), values(i), 1, 'Hazen');
        end
end

