function [pInds, k_adjusted] = uperms(X, k, ulim)
%uperms: unique permutations of an input vector or rows of an input matrix
% Usage:  nPerms              = uperms(X)
%        [nPerms pInds]       = uperms(X, k)
%        [nPerms pInds Perms] = uperms(X, k)
%
% Determines number of unique permutations (nPerms) for vector or matrix X.
% Optionally, all permutations' indices (pInds) are returned. If requested,
% permutations of the original input (Perms) are also returned.
%
% If k < nPerms, a random (but still unique) subset of k of permutations is
% returned. The original/identity permutation will be the first of these.
%
% Row or column vector X results in Perms being a [k length(X)] array,
% consistent with MATLAB's built-in perms. pInds is also [k length(X)].
%
% Matrix X results in output Perms being a [size(X, 1) size(X, 2) k]
% three-dimensional array (this is inconsistent with the vector case above,
% but is more helpful if one wants to easily access the permuted matrices).
% pInds is a [k size(X, 1)] array for matrix X.
%
% Note that permutations are not guaranteed in any particular order, even
% if all nPerms of them are requested, though the identity is always first.
%
% Other functions can be much faster in the special cases where they apply,
% as shown in the second block of examples below, which uses perms_m.
%
% Examples:
%  uperms(1:7),       factorial(7)        % verify counts in simple cases,
%  uperms('aaabbbb'), nchoosek(7, 3)      % or equivalently nchoosek(7, 4).
%  [n pInds Perms] = uperms('aaabbbb', 5) % 5 of the 35 unique permutations
%  [n pInds Perms] = uperms(eye(3))       % all 6 3x3 permutation matrices
%
%  % A comparison of timings in a special case (i.e. all elements unique)
%  tic; [nPerms P1] = uperms(1:20, 5000); T1 = toc
%  tic; N = factorial(20); S = sample_no_repl(N, 5000);
%  P2 = zeros(5000, 20);
%  for n = 1:5000, P2(n, :) = perms_m(20, S(n)); end
%  T2 = toc % quicker (note P1 and P2 are not the same random subsets!)
%  % For me, on one run, T1 was 7.8 seconds, T2 was 1.3 seconds.
%
%  % A more complicated example, related to statistical permutation testing
%  X = kron(eye(3), ones(4, 1));  % or similar statistical design matrix
%  [nPerms pInds Xs] = uperms(X, 5000); % unique random permutations of X
%  % Verify correctness (in this case)
%  G = nan(12,5000); for n = 1:5000; G(:, n) = Xs(:,:,n)*(1:3)'; end
%  size(unique(G', 'rows'), 1)    % 5000 as requested.
%
% See also: randperm, perms, perms_m, signs_m, nchoosek_m, sample_no_repl
% and http://www.fmrib.ox.ac.uk/fsl/randomise/index.html#theory
% Copyright 2010 Ged Ridgway
% http://www.mathworks.com/matlabcentral/fileexchange/authors/27434
%% Usage
if isempty(X) || isscalar(X) || ~ismatrix(X)
    error('Input X must be a vector or matrix')
end
if nargin > 1 && nargout < 2
    error('Set number of permutations requested, but no return argument!')
end
if islogical(X), X = find(X); end

if ~exist("ulim","var") || isempty(ulim)
    if isvector(X)
        ulim = size(X,1);
    else
        ulim = size(X,2);
    end
end
k_adjusted = k;
%% Count number of repetitions of each unique row, and get representative x
if isvector(X)
    [~, ~, x] = unique(X(:)); % x codes unique elements with integers
else
    [~, ~, x] = unique(X, 'rows'); % x codes unique rows with integers
end
nPerms = nchoosek(numel(x), ulim);

%% Computation of permutations
% Basics
n = length(x);
if nPerms < k 
    k = nPerms; 
    k_adjusted = k;
end
if n > uint32(inf), error('Sorry, data is too big!'), end % would be v.slow
if nargin < 2 || k > nPerms
    k = nPerms; % default to computing all unique permutations
end
pInds = false(k, n);              
% Add permutations that are unique
u = 1; % to start with
fprintf('Creating unique perms:    \n\n')
while u <= k
    fprintf('\b\b\b\b\b\b%6.0f',u);
    % b=sprintf('%s', repmat('\b',numel(a),1));
    pInd = sort(uint32(randperm(n, ulim)));
    % (Note: better data structures could make the next line much faster!)
    tpInd = false(1, n); tpInd(pInd) = true;
    out = FastUniqueRows([pInds(1:u, :); tpInd]);
    if size(out, 1) > u % implying x(pInd) was new to Perms
        pInds(u, pInd) = true ;
        u = u + 1;
    end
end