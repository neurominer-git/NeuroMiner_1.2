function [rs, ds] = nk_GetTestPerf_DECTRE(~, tXtest, ~, md, ~, ~)

[rs, ds] = md.predict(tXtest); 
ds = ds(:,2)-0.5;


