function [sY, sX, sYX, IN] = nk_PLS(Y, X, IN)
% Computes PLS model from experimental predictors X and predictands Y
% _________________________________________________________________________
% (c) Nikolaos Koutsouleris, 08/2023
if ~isfield(IN,'trained'),IN.trained=0;end
if ~isfield(IN,'algostr'),IN.algostr='pls';end
% Get behavioural data
if exist('X','var') && ~isempty(X)
    Xd = X;
else
    Xd = IN.DR.PLS.V;
end

if ~IN.trained

    % Mean centering of Y and Xd
    switch IN.algostr
        case 'pls'
            IN.mpp.mY.method = 'mean-centering';
            IN.mpp.mX.method = 'mean-centering';
        case 'spls'
            IN.mpp.mY.method = 'standardization using mean';
            IN.mpp.mX.method = 'standardization using mean';
    end
    [mY, IN.mpp.mY] = nk_PerfStandardizeObj(Y,IN.mpp.mY);
    [mXd, IN.mpp.mX] = nk_PerfStandardizeObj(Xd,IN.mpp.mX);
    
    % If mXd contains NaNs, we have to implement an imputation step here:
    % Imputation
    if sum(isnan(mXd(:)))
        mXd = SeqkNN(mXd, 5);
        IN.tX = mXd;
    end

    switch IN.algostr
        case 'pls'
            % Build covariance matrix S
            S=mXd'*mY;
            % Decompose S into LVpairs
            [IN.mpp.u, IN.mpp.s, IN.mpp.v] = svd(S',0);
            IN.mpp.d = diag(IN.mpp.s)';
        case 'spls'
            % Perform sparse PLS
            nD = size(mXd,2);
            IN.mpp.u = zeros(size(mY,2),nD);
            IN.mpp.v = zeros(nD,nD);
            IN.mpp.C = [];
            for i=1:nD
                % Compute SPLS level i
                [IN.mpp.u(:,i), IN.mpp.v(:,i), IN.mpp.C] = spls(mY, mXd, IN.cu, IN.cv, IN.mpp.C);
                % Covariance matrix deflation
                IN.mpp.C = IN.mpp.C - (IN.mpp.C * IN.mpp.v(:,i)') * IN.mpp.v(:,i)';
            end
    end
    sY = mY*IN.mpp.u;
    sX = mXd*IN.mpp.v;
    sYX = mXd * (IN.mpp.u * IN.mpp.v')' + IN.mpp.mY.meanY;
    IN.trained = true;
else
    mY = nk_PerfStandardizeObj(Y,IN.mpp.mY);
    sY = mY*IN.mpp.u; sX = []; sYX = [];
    if ~isempty(Xd)
        mXd = nk_PerfStandardizeObj(Xd,IN.mpp.mX);
        % If mXd contains NaNs, we have to implement an imputation step here:
        % Imputation
        if sum(isnan(mXd(:)))
            mXd = SeqkNN(mXd, 5, IN.tX);
        end
        sX = mXd*IN.mpp.v;
        sYX = mXd * (IN.mpp.u * IN.mpp.v')' + IN.mpp.mY.meanY;
    end
end
