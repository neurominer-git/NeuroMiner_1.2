% Decile-based plot of probabilities and empirical frequencies of labels
% If the two sit on top of each other, that indicates the probabilities are
% well-calibrated
% p is a real-valued vector of probability scores
% y is a binary vector of labels
function plotCalibrationDeciles(p, y, deciles)

    [v,i] = sort(p);
    p = p(i);
    y = y(i);

    % How many deciles?
    if nargin < 3, deciles = 100; end;

    N = floor(length(y)/deciles) * deciles;
    yReshaped = reshape(y(1:N)',[N/deciles deciles]);
    pReshaped = reshape(p(1:N),[N/deciles deciles]);
    figure;
    plot(mean(pReshaped), mean(yReshaped))
