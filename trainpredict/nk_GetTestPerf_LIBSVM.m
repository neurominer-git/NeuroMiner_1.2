function [rs, ds] = nk_GetTestPerf_LIBSVM(~, tXtest, Ytest, md, ~, ~)
global SVM LIBSVMPREDICT

% Check whether a posthoc calibration model needs to be applied to test estimates
if isfield(md,'BBQ') && ~SVM.LIBSVM.Optimization.b
    [~, ~, PTE] = feval( LIBSVMPREDICT, Ytest, tXtest, md.md, sprintf(' -b %g',SVM.LIBSVM.Optimization.b)); 
    PTE = exp(PTE)./(1+exp(PTE)); predict_test_calib = predictBBQ(md.BBQ,PTE,0); 
    rs = sign(predict_test_calib-.5); ds = predict_test_calib;
else
    [err_test, ~, predict_test] = feval( LIBSVMPREDICT, Ytest, tXtest, md, sprintf(' -b %g',SVM.LIBSVM.Optimization.b)); 
    rs = err_test; ds = predict_test(:,1);
end

