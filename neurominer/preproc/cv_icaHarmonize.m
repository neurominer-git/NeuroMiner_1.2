function [ Y, IN ] = cv_icaHarmonize(Y, mode, IN)

if mode == 1

    if isfield(IN,'subgroup') && ~isempty(IN.subgroup)
        Ytemp = Y; 
        Y = Y(IN.subgroup,:);
        G = IN.G(IN.subgroup,:);
    else
        G = IN.G; 
    end

    [icasig, A, W] = fastica(Y');
    % icasig = projection of original data on the ICs
    % A = ICs are represented in the columns
    % W = separating matrix W can be used to project new data on ICs (Yn =
    % W*Xn)

    % Correlate each independent component with the covariate
    correlation = zeros(1, size(icasig, 1));
    pvals = zeros(1, size(icasig, 1));
    for i = 1:size(icasig, 1)
        [correlation(i), pvals(i)] = corr(icasig(i, :).', G);
    end

    % Identify the independent components significantly associated with the covariate
    if ~isfield(IN, "pthres")
        IN.pthres = 0.05;
    end
    significantComponents = find(abs(pvals) < IN.pthres);


    if isfield(IN,'subgroup') && ~isempty(IN.subgroup)
        % Remove the significant components from the original data
        Ysub = Y - (W(significantComponents, :)' * icasig(significantComponents, :))';

        % project left out data on ICs
        Ysubinv = Ytemp(~IN.subgroup,:);
        Yn = W * Ysubinv';
        % Remove the significant components from the second dataset
        Ysubinv = Ysubinv - (W(significantComponents, :)' * Yn(significantComponents, :))';

        Ytemp(IN.subgroup,:) = Ysub;
        Ytemp(~IN.subgroup,:) = Ysubinv;
        Y = Ytemp;
    else
        Y = Y - (W(significantComponents, :)' * icasig(significantComponents, :))';
    end

    IN.W = W;
    IN.icasig = icasig;
    IN.A = A;
    IN.sigICs = significantComponents;
    IN.pvals = pvals;
else
    % project new data on ICs
    Yn = IN.W * Y';

    % Remove the significant components from the second dataset
    Y = Y - (IN.W(IN.sigICs, :)' * Yn(IN.sigICs, :))';

end