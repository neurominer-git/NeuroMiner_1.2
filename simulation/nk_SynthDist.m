% Load the real training data and labels
load('real_data.mat');

% Set the number of synthetic observations to generate
numSyntheticObservations = 1000;

% Initialize variables to store the synthetic data and labels
syntheticData = zeros(numSyntheticObservations, size(realData, 2));
syntheticLabels = zeros(numSyntheticObservations, 1);

% Loop through each real observation
for i = 1:numSyntheticObservations
    % Select a random real observation to use as the basis for the synthetic observation
    realObservationIndex = randi([1, size(realData, 1)]);
    realObservation = realData(realObservationIndex, :);
    realLabel = realLabels(realObservationIndex, :);

    % Generate a synthetic observation by adding random noise to the selected real observation
    sigma = 1;
    noise = random('Normal', 0, sigma, 1, size(realData, 2));
    syntheticObservation = realObservation + noise;

    % Compute the distances between the synthetic observation and all real observations
    distances = sqrt(sum((realData - repmat(syntheticObservation, size(realData, 1), 1)).^2, 2));

    % Select the closest real observation to the synthetic observation as the label for the synthetic observation
    [~, closestRealObservationIndex] = min(distances);
    syntheticLabel = realLabels(closestRealObservationIndex, :);

    % Store the synthetic observation and label in the appropriate arrays
    syntheticData(i, :) = syntheticObservation;
    syntheticLabels(i, :) = syntheticLabel;
end

% Combine the real and synthetic data and labels
data = [realData; syntheticData];
labels = [realLabels; syntheticLabels];

% =========================================================================
% Load the real training data and labels
load('real_data.mat');

% Set the number of synthetic observations to generate
numSyntheticObservations = 1000;

% Set the number of closest neighbors to use for interpolation
k = 5;

% Set the distance measure to use (options: 'Euclidean', 'Manhattan', 'Cosine', 'Mahalanobis', 'Wasserstein')
distanceMeasure = 'Euclidean';

% Initialize variables to store the synthetic data and labels
syntheticData = zeros(numSyntheticObservations, size(realData, 2));
syntheticLabels = zeros(numSyntheticObservations, 1);

if strcmp(distanceMeasure,'Mahalanobis')
    % Compute the covariance matrix of the real data (needed for the Mahalanobis distance)
    covarianceMatrix = cov(realData);
end

% Loop through each real observation
for i = 1:numSyntheticObservations
    % Select a random real observation to use as the basis for the synthetic observation
    realObservationIndex = randi([1, size(realData, 1)]);
    realObservation = realData(realObservationIndex, :);
    realLabel = realLabels(realObservationIndex, :);

    % Compute the distances between the selected real observation and all other real observations
    switch distanceMeasure
        case 'Euclidean'
            distances = sqrt(sum((realData - repmat(realObservation, size(realData, 1), 1)).^2, 2));
        case 'Manhattan'
            distances = sum(abs(realData - repmat(realObservation, size(realData, 1), 1)), 2);
        case 'Cosine'
            distances = 1 - sum((realData .* repmat(realObservation, size(realData, 1), 1)), 2) ./ (sqrt(sum(realData.^2, 2)) .* sqrt(sum(realObservation.^2)));
        case 'Mahalanobis'
            distances = sqrt(mahal(realData, realObservation, covarianceMatrix));
        case 'Wasserstein'
            distances = wasserstein(realData, realObservation);
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
end

% Combine the real and synthetic data and labels
data = [realData; syntheticData];
labels = [realLabels; syntheticLabels];

% Train the kNN model using the combined data
knn = fitcknn(data, labels, 'NumNeighbors', 5);

% Compute the Wasserstein distance between two vectors
function distance = wasserstein_gpt3(x, y)
    % Compute the distance matrix between the two vectors
    distanceMatrix = abs(repmat(x, size(y, 1), 1) - repmat(y, size(x, 1), 1));

    % Compute the Wasserstein distance as the sum of the elements of the distance matrix divided by the number of elements
    distance = sum(distanceMatrix(:)) / numel(distanceMatrix);
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
u_samples_sorted = sort(u_samples(:));
v_samples_sorted = sort(v_samples(:));
if p == 1
    
    all_samples = unique([u_samples_sorted; v_samples_sorted], 'sorted');
    
    u_cdf = find_interval(u_samples_sorted, all_samples(1:end-1)) ...
        / numel(u_samples);
    v_cdf = find_interval(v_samples_sorted, all_samples(1:end-1)) ...
        / numel(v_samples);
    
    wsd = sum(abs(u_cdf - v_cdf) .* diff(all_samples));
    
elseif p == 2
    
    u_N = numel(u_samples);
    v_N = numel(v_samples);    
    all_prob = unique([(0:u_N) / u_N, (0:v_N) / v_N], 'sorted').';
    
    u_icdf = u_samples_sorted(fix(all_prob(1:end-1) * u_N) + 1);
    v_icdf = v_samples_sorted(fix(all_prob(1:end-1) * v_N) + 1);
    
    wsd = sqrt(sum((u_icdf-v_icdf).^2 .* diff(all_prob)));
    
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
end

