function [data, labels, covars] = nk_SynthDistkNN(realData, realLabels, realCovars, IN)

if ~exist("IN","var") || isempty(IN)
    distanceMeasure = 'Euclidean'; 
    k = 5;
    numSyntheticObservations = 100;
    std_flag = 2;
else
    distanceMeasure = IN.distanceMeasure;
    k = IN.k;
    numSyntheticObservations = IN.numSyntheticObservations;
    std_flag = IN.standardize_data;
end

% Initialize variables to store the synthetic data and labels
syntheticData = zeros(numSyntheticObservations, size(realData, 2));
syntheticLabels = zeros(numSyntheticObservations, 1);
if ~isempty(realCovars)
    syntheticCovars = zeros(numSyntheticObservations, size(realCovars,2));
else
    syntheticCovars = [];
end
fprintf('\nGenerating synthetic data using kNN-based random interpolation.');
fprintf('\n\tDistance Measure: %s', distanceMeasure);
fprintf('\n\tk: %g',k);
fprintf('\n\tNo. of synthetic observations: %g\n', numSyntheticObservations);

% Do we need to standardize the data?
if std_flag == 1
    fprintf('\nStandardizing data')
    std_realData = nk_PerfStandardizeObj(realData);
else
    std_realData = realData;
end

% Check whether there are NaNs in the data and impute them
if sum(isnan(std_realData(:)))
    fprintf('\nImputing missing values in data');
    IN_impute.k = k;
    IN.impute.method = 'SeqkNN'; 
    std_realData = nk_PerfImputeObj(std_realData, IN_impute);
end

if strcmp(distanceMeasure,'Mahalanobis')
    % Compute the covariance matrix of the real data (needed for the Mahalanobis distance)
    covarianceMatrix = cov(realData);
end


% Initialize a logical array to keep track of used real observations
usedRealObservations = false(size(realData, 1), 1);

% Loop through each real observation
for i = 1:numSyntheticObservations

    % Select a random unused real observation to use as the basis for the synthetic observation
    unusedRealObservationIndices = find(~usedRealObservations);
    if isempty(unusedRealObservationIndices)
        warning('No unused real observations left to create synthetic data.');
        break
    end
    % Select a random real observation to use as the basis for the synthetic observation
    randomIndex = randi(length(unusedRealObservationIndices));
    realObservationIndex = unusedRealObservationIndices(randomIndex);
    usedRealObservations(realObservationIndex) = true;
    std_realObservation = std_realData(realObservationIndex, :);
    
    fprintf('.')
    % Compute the distances between the selected real observation and all other real observations
    switch distanceMeasure
        case 'Euclidean'
            distances = sqrt(sum((std_realData - repmat(std_realObservation, size(realData, 1), 1)).^2, 2));
        case 'Manhattan'
            distances = sum(abs(std_realData - repmat(std_realObservation, size(realData, 1), 1)), 2);
        case 'Cosine'
            distances = 1 - sum((std_realData .* repmat(std_realObservation, size(realData, 1), 1)), 2) ./ (sqrt(sum(std_realData.^2, 2)) .* sqrt(sum(std_realObservation.^2)));
        case 'Mahalanobis'
            distances = sqrt(mahal(std_realData, std_realObservation, covarianceMatrix));
        case 'Wasserstein'
            distances = wasserstein(realData, std_realObservation);
    end

    % Select the k closest real observations to the selected real observation
    [~, closestRealObservationIndices] = mink(distances, k);
    closestRealObservations = realData(closestRealObservationIndices, :);
    closestRealLabels = realLabels(closestRealObservationIndices, :);
    
    % Generate a synthetic observation by interpolating between the k closest real observations
    syntheticObservation = mean(closestRealObservations, 1);

    % Generate a synthetic label by using the mode of the k closest real labels
    syntheticLabel = mode(closestRealLabels, 1);
    
    % Store the synthetic observation and label in the appropriate arrays
    syntheticData(i, :) = syntheticObservation;
    syntheticLabels(i, :) = syntheticLabel;
    
    if ~isempty(realCovars)
        closestRealCovars = realCovars(closestRealObservationIndices, :);
        syntheticCovar = mean(closestRealCovars, 1);
        syntheticCovars(i,:) = syntheticCovar;
    end
end
fprintf('\nDone!\n')
% Combine the real and synthetic data and labels
data = [realData; syntheticData];
labels = [realLabels; syntheticLabels];
if ~isempty(realCovars)
    covars = [realCovars; syntheticCovars];
end

function wsd = wasserstein(u_samples, v_samples, p)
% WS_DISTANCE 1- and 2- Wasserstein distance between two discrete 
% probability measures 
%   
%   wsd = WS_DISTANCE(u_samples, v_samples) returns the 1-Wasserstein 
%   distance between the discrete probability measures u and v 
%   corresponding to the sample vectors u_samples and v_samples
%
%   wsd = WS_DISTANCE(u_samples, v_samples, p) returns the p-Wasserstein 
%   distance between the discrete probability measures u and v
%   corresponding to the sample vectors u_samples and v_samples. 
%   p must be 1 or 2.
%
% from https://github.com/nklb/wasserstein-distance
if ~exist('p', 'var')
    p = 1;
end

uN = size(u_samples,1);
wsd=zeros(uN,1);
v_samples_sorted = sort(v_samples(:));
u_samples=u_samples';
for i=1:uN

    u_samples_sorted = sort(u_samples(:,i));

    if p == 1
        
        all_samples = unique([u_samples_sorted; v_samples_sorted], 'sorted');
        
        u_cdf = find_interval(u_samples_sorted, all_samples(1:end-1)) ...
            / numel(u_samples);
        v_cdf = find_interval(v_samples_sorted, all_samples(1:end-1)) ...
            / numel(v_samples);
        
        wsd(i) = sum(abs(u_cdf - v_cdf) .* diff(all_samples));
        
    elseif p == 2
        
        u_N = numel(u_samples);
        v_N = numel(v_samples);    
        all_prob = unique([(0:u_N) / u_N, (0:v_N) / v_N], 'sorted').';
        
        u_icdf = u_samples_sorted(fix(all_prob(1:end-1) * u_N) + 1);
        v_icdf = v_samples_sorted(fix(all_prob(1:end-1) * v_N) + 1);
        
        wsd(i) = sqrt(sum((u_icdf-v_icdf).^2 .* diff(all_prob)));
        
    else
        
        error('Only p=1 or p=2 allowed.')
        
    end
end

function idx = find_interval(bounds, vals)
% Given the two sorted arrays bounds and vals, the function 
% idx = FIND_INTERVAL(bounds, vals) identifies for each vals(i) the index 
% idx(i) s.t. bounds(idx(i)) <= vals(i) < bounds(idx(i) + 1).
m = 0;
bounds = [bounds(:); inf];
idx = zeros(numel(vals), 1);
for i = 1:numel(vals)
    while bounds(m+1) <= vals(i)
        m = m + 1;
    end
    idx(i) = m;
end


