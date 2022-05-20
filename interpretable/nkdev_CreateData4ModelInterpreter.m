function [P, N, I] = nkdev_CreateData4ModelInterpreter(Tr, Ts, nperms, frac, mapthresh)

% Determine extremes of the distribution
[~, n] = size(Tr,2);
upper = percentile(Tr, 95);
lower = percentile(Tr, 5);
P = repmat(Ts, nperms, 1);
N = repmat(Ts, nperms, 1);

% If map is provided determine subspace for modification
if ~isempty(mapthresh)
    mapidx = return_imgind(mapthresh.typ, mapthresh.thresh, mapthresh.map);
else
    mapidx = 1:n;
end

n = numel(mapidx);
nfrac = ceil(n*frac);
I = zeros(nperms, nfrac);

% Create modified instances of case
for i=1:nperms
    idx = randperm(nfrac, n);
    P(i, idx) = upper(idx); 
    N(i, idx) = lower(idx);
    I(i, :) = idx;
end

function imgind = return_imgind(typthresh, thresh, img)

if length(thresh) > 1
    switch typthresh
        case 1
            imgind = (img < thresh(1) | img > thresh(2)); 
        case 2
            imgind = (img <= thresh(1) | img >= thresh(2)); 
        case 3
            imgind = (img > thresh(1) | img < thresh(2)); 
        case 4
            imgind = (img >= thresh(1) | img <= thresh(2)); 
    end
else
    switch typthresh
        case 1
            imgind = img < thresh; 
        case 2
            imgind = img <= thresh;
        case 3
            imgind = img > thresh;
        case 4
            imgind = img >= thresh;
        case 5
            imgind = img == thresh;
    end
end
imgind = find(imgind);