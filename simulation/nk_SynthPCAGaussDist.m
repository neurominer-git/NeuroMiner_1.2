% Load the real training data and labels
load('real_data.mat');

% Set the number of synthetic observations to generate
numSyntheticObservations = 1000;

% Perform principal component analysis on the real data
[coeff, score, latent] = pca(realData);

% Generate synthetic data and labels by sampling random values from a Gaussian distribution
mu = mean([score, realLabels], 1);
sigma = std([score, realLabels], 1);
syntheticDataAndLabels = random('Normal', repmat(mu, numSyntheticObservations, 1), repmat(sigma, numSyntheticObservations, 1));

% Project the synthetic data onto the principal components of the real data
syntheticData = syntheticDataAndLabels(:, 1:size(realData, 2)) * coeff';

% Extract the synthetic labels
syntheticLabels = syntheticDataAndLabels(:, size(realData, 2)+1:end);

% Combine the real and synthetic data and labels
data = [realData; syntheticData];
labels = [realLabels; syntheticLabels];

% Train the SVM using the combined data
svm = fitcsvm(data, labels, 'KernelFunction', 'linear', 'BoxConstraint', 1);

