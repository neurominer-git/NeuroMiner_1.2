function IN = nk_PerfANOVAObjNew(Y, IN)
global VERBOSE
[~,n] = size(Y);
%IN.p = zeros(n,1);
%IN.F = zeros(n,1);
%IN.R2 = zeros(n,1);
if VERBOSE, fprintf('\t running ANOVA on %g variables ',n); end

IN.beta = pinv(IN.X)*Y; 
% here we use a pseudoinverse:
% X is rank deficient, i.e. regressors are not independent, since any linear combination of
% 3 columns can give us the 4th one, thus X'X is also rank deficient and singular ie inv(X'*X)
% doesn't exist -- there is no matrix A such as A(X'X) = I - however pinv gives a unique
% solution that minimizes the square distance to the data.
% see
% http://mathworld.wolfram.com/Moore-PenroseMatrixInverse.html
% http://en.wikipedia.org/wiki/Moore_Penrose_pseudoinverse

IN.Yhat = IN.X*IN.beta;
IN.df = rank(IN.X)-1;

IN.Res = Y - IN.Yhat;
IN.SStotal = vecnorm(Y - nm_nanmean(Y)).^2;
IN.SSeffect = vecnorm(IN.Yhat - nm_nanmean(IN.Yhat)).^2;
IN.SSerror  = vecnorm(IN.Res - nm_nanmean(IN.Res)).^2;
IN.dferror = size(Y,1) - IN.df - 1;
IN.R2 = IN.SSeffect ./ IN.SStotal;
IN.F = (IN.SSeffect / IN.df) ./ (IN.SSerror / IN.dferror);
IN.p = 1 - fcdf(IN.F,IN.df,IN.dferror);

