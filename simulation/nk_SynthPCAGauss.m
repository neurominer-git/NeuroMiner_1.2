function [data, labels, covars] = nk_SynthPCAGauss(realData, realLabels, realCovars, IN)

if ~exist("IN","var") || isempty(IN)
    numSyntheticObservations = 100;
else
    numSyntheticObservations = IN.numSyntheticObservations ;
end

fprintf('\nGenerating synthetic data using Gaussian random interpolation.');
fprintf('\n\tNo. of synthetic observations: %g\n', numSyntheticObservations);

dRealData = size(realData,2);
if ~isempty(realCovars)
    dRealCovars = size(realCovars,2);
    realData= [realData realCovars realLabels];
    uXC = nk_CountUniques(realCovars);
else
    realData = [realData realLabels];
    dRealCovars = 0;
end

uXL = nk_CountUniques(realLabels);

% Check whether there are NaNs in the data and impute them
if sum(isnan(realData(:)))
    fprintf('\nImputing missing values in data');
    IN_impute.k = 5;
    IN_impute.method = 'SeqkNN'; 
    IN_impute.X = realData;
    realData = nk_PerfImputeObj(realData, IN_impute);
end

% Perform principal component analysis on the combined data
IN.DR.DRsoft = 1; 
IN.DR.RedMode = 'PCA';
IN.DR.PercMode = 3; 
[score, model] = nk_PerfRedObj(realData, IN);

% Generate synthetic data and labels by sampling random values from a Gaussian distribution
mu = mean(score, 1);
sigma = std(score, 1);
syntheticDataCovarsLabels = random('Normal', repmat(mu, numSyntheticObservations, 1), repmat(sigma, numSyntheticObservations, 1));

% Project the synthetic data onto the principal components of the combined data
tsyntheticDataCovarsLabels = syntheticDataCovarsLabels * model.mpp.vec' + model.mpp.sampleMean;
syntheticDataCovarsLabels = zeros(numSyntheticObservations,numel(model.indNonRem));
syntheticDataCovarsLabels(:,model.indNonRem) = tsyntheticDataCovarsLabels ;
clear tsyntheticDataCovarsLabels

% Extract synthetic data
syntheticData = syntheticDataCovarsLabels(:,1:dRealData);

% Extract covariates
if ~isempty(realCovars)
    syntheticCovars = syntheticDataCovarsLabels(:,dRealData+1:dRealData+dRealCovars);
    rIdx = find(uXC.U<10);
    if ~isempty(rIdx)
        for i=1:numel(rIdx)
            syntheticCovars(:,rIdx(i)) = interp1(uXC.UX{rIdx(i)},uXC.UX{rIdx(i)},syntheticCovars(:,rIdx(i)),'nearest','extrap');
        end
    end
end

% Extract the synthetic labels
syntheticLabels = syntheticDataCovarsLabels(:, dRealData+dRealCovars+1:end);
rIdx = find(uXL.U<10);
if ~isempty(rIdx)
    for i=1:numel(rIdx)
        syntheticLabels(:,rIdx(i)) = interp1(uXL.UX{rIdx(i)},uXL.UX{rIdx(i)},syntheticLabels(:,rIdx(i)),'nearest','extrap');
    end
end

fprintf('\nDone!\n')
% Combine the real and synthetic data and labels
data = syntheticData;
labels = syntheticLabels;
if ~isempty(realCovars)
    covars = syntheticCovars;
end
