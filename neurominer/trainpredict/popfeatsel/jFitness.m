% This fitness function for feature selection
function cost = jFitnessFunction1(r , X, Ps)
if sum(X == 1) == 0
  cost = 1;
else
  cost = jwrapper(r.Y(:, X), r.YL, r.T(:, X), r.L, Ps );
end

function cost = jwrapper(xtrain, ytrain, xvalid, yvalid, Ps)
global TRAINFUNC

[~, FullModel] = TRAINFUNC(xtrain, ytrain, 1, Ps);    
cost = nk_GetTestPerf(xvalid, yvalid, [], FullModel, xtrain);


