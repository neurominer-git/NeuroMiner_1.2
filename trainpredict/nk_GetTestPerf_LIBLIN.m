function [rs, ds] = nk_GetTestPerf_LIBLIN(~, tXtest, Ytest, md, ~, ~)
global SVM

% Check whether a posthoc calibration model needs to be applied to test estimates
if isfield(md,'BBQ') && ~SVM.LIBLIN.b 
    [~, ~, PTE] = predict_liblin22(Ytest, sparse(tXtest), md.md, sprintf(' -b %g -q',SVM.LIBLIN.b)); 
    PTE = exp(PTE(:,1))./(1+exp(PTE(:,1))); predict_test_calib = predictBBQ(md.BBQ,PTE,0);  
    rs = sign(predict_test_calib-.5); ds = predict_test_calib;
else
    [err_test, ~, predict_test] = predict_liblin22(Ytest, sparse(tXtest), md, sprintf(' -b %g -q',SVM.LIBLIN.b )); 
    rs = err_test; ds = predict_test(:,1);
end