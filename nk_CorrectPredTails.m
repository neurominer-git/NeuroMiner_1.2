function [ sY, Y, IN ] = nk_CorrectPredTails( Y, X, IN )
% =========================================================================
% function [ sY, Y, IN ] = nk_CorrectPredTails( Y, X, IN )
% =========================================================================
% Regression models usually produce higher prediction errors at the tails
% of the label distribution. This functon corrects for this effect by
% (1) computing the slope between labels and prediction errors using a
% reference dataset and (2) applying the correction parameters to the
% predictions of a target sample.
%
% IN could be defined as follows, e.g.:
% IN.TrPred = NM.analysis{1}.GDdims{1}.Regr.mean_predictions
% IN.TrObs = NM.label;
%
% Y could be :
% NM.analysis{1}.OOCV{1}.RegrResults{1}.Group{1}.MeanCV2PredictedValues
% _________________________________________________________________________
% (c) Nikolaos Koutsouleris, 02/2021

flag = false;
if exist('IN','var') && ~isempty(IN) 
    if ~isfield(IN,'beta') || ~isfield(IN,'p')
        if ~isfield(IN,'TrObs') || ~isfield(IN,'TrPred')
            error('Please provide the reference sample''s labels and predictions so that I can compute correction parameters!')
        else
            if isempty(Y)
                % LOO mode, if no target sample has been provided with Y
                flag = true;
                nTr = numel(IN.TrPred);
                Y = IN.TrPred;
                sY = zeros(nTr,1);
                IN.beta = zeros(2,nTr);
                IN.p =zeros(nTr,2);
                fprintf('\nLOO mode')
                for i=1:nTr
                    fprintf('.')
                    I_train = true(nTr,1); I_train(i)=false;
                    [~, IN.beta(:,i), IN.p(i,:)] = nk_DetrendPredictions2([],[], Y(I_train), IN.TrObs(I_train));
                    sY(i) = nk_DetrendPredictions2(IN.beta(:,i),IN.p(i,:), Y(i));
                end
                fprintf('\n')
            else
                % Reference modelling mode
                [~, IN.beta, IN.p] = nk_DetrendPredictions2([],[], IN.TrPred, IN.TrObs);
            end
        end
    end
else
    error('Please provide a valid IN structure for the function!')
end
if ~exist("X", "var"), X=[]; end
if ~flag, sY = nk_DetrendPredictions2(IN.beta, IN.p, Y, X); end
