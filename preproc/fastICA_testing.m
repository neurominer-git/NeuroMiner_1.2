%% Prep: Generate example data
rng(123); % Set the random seed
n = 1000;
x = randn(n, 3);
A = [2, 1, -1; 1, 2, 1; 1, -1, 2];
s = x * A'; % Mix the signals

% in the Matlab implementation of fastica(), whitening of the data is
% internally performed, so we do not need to do it before hand 
% Perform ICA using FastICA
[ica.S, ica.A, ica.W] = fastica(s', 'numOfIC', 3); % default for numOfIC = min(size(s))

%% 1) Compute the original data from the independent components
reconstructed_data = ica.S' * ica.A;

% View the results
figure;
subplot(2,1,1); plot(x); title('Original Data');
subplot(2,1,2); plot(reconstructed_data); title('Reconstructed Data from ICA Output');


% Note that the reconstructed data may differ from the original data due to the 
% fact that ICA is an unsupervised method that does not necessarily find the exact
% same mixing matrix that was used to generate the mixed signals.
% 
% ICA is an unsupervised method that relies on statistical independence to separate 
% mixed signals into their independent components. The algorithm searches for a 
% set of uncorrelated signals that have non-Gaussian distributions, which are more 
% likely to be statistically independent. Because ICA is an unsupervised method, 
% it does not know the true mixing matrix that was used to create the mixed signals. 
% Instead, it estimates the mixing matrix based on the statistical properties of 
% the observed mixed signals. This means that the estimated mixing matrix may not 
% be exactly the same as the true mixing matrix, especially if the observed signals 
% are noisy or if there are other sources of variation that the algorithm cannot 
% account for. Furthermore, even if the estimated mixing matrix is close to the 
% true mixing matrix, there may still be some residual error in the reconstructed 
% data due to the fact that the algorithm is based on statistical estimation rather 
% than exact computation. Therefore, the reconstructed data from the ICA output may 
% not be exactly the same as the original data, but it should be a reasonable 
% approximation that captures the underlying independent components that make up 
% the original signals.

%% 2) Generate new data and project onto independent components
newdata = randn(n, 3); % New data
newdata_centered = newdata - mean(x); % Center the new data
newdata_whitened = newdata_centered * ica.W'; % Whiten the new data
newdata_proj = newdata_whitened * ica.A; % Project the new data onto the independent components using the same weights as the training data

% Plot the results
figure;
subplot(2,1,1);
plot(newdata);
title('New Data');
subplot(2,1,2);
plot(newdata_proj);
title('New Data Projected onto Found Components');


%% 3) fastICA to reduce dimensionality and project testset on ICs extracted in trainset
% Generate example data
rng(123); % Set random seed
n = 1000;
x_train = randn(n, 5); % Original training data with 5 features
x_test = randn(500, 5); % Original test data with 5 features
A = [2 1 -1 0 0; 1 2 1 0 0; 1 -1 2 0 0; 0 0 0 2 1; 0 0 0 1 2]; % Mixing matrix
s_train = x_train * A'; % Mix the training signals
s_test = x_test * A'; % Mix the test signals

% Perform ICA using FastICA
[icasig_train, A, W] = fastica(s_train', 'approach', 'symm', 'numOfIC', 3);

% Use the independent components as features for classification
X_train = icasig_train'; % Transpose icasig_train to have features as rows and samples as columns
Y_train = [ones(n/2,1); 2*ones(n/2,1)]; % Create binary labels for a two-class problem

% Train a classifier (example: linear SVM)
mdl = fitcsvm(X_train, Y_train);

% Extract independent components from test data
icasig_test = W * s_test';
X_test = icasig_test'; % Transpose icasig_test to have features as rows and samples as columns

% Make predictions on test data using the trained classifier
Y_test = [ones(500/2,1); 2*ones(500/2,1)]; % Create binary labels for a two-class problem
Y_pred = predict(mdl, X_test);

% Evaluate the accuracy of the classifier on test data
accuracy = sum(Y_pred == Y_test) / length(Y_test);
fprintf('Accuracy: %.2f%%\n', 100*accuracy);