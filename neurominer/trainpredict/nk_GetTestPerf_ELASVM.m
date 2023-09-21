function [rs, ds] = nk_GetTestPerf_ELASVM(~, tXtest, ~, model, ~, ~)

z = tXtest *model.beta; ds = 1 ./ (1 + exp(-z)); ds = ds - 0.5;
rs = sign(ds); 