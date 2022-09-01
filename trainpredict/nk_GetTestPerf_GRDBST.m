function [rs, ds] = nk_GetTestPerf_GRDBST(~, tXtest, ~, md, ~, ~)
global MODEFL MODELDIR

%ds = SQBMatrixPredict( md, single(tXtest));
switch MODEFL
    case 'classification'
        %[rs, votes] = classRF_predict(tXtest,md);
        results_file = pyrunfile('py_classGRDBST_predict.py', ...
            'results_file' , model_name = md, test_feat =tXtest, ...
            rootdir = MODELDIR);
        results = load(char(results_file));
        rs = results.predictions;
        votes = results.probabilities;
        ds = votes(:,2)./sum(votes,2);
    case 'regression'
        %rs = regRF_predict(tXtest,md); ds=rs;
        results_file = pyrunfile('py_regGRDBST_predict.py', ...
            'results_file', model_name = md, teast_feat = tXtest, ...
            rootdir = MODELDIR);
        results = load(char(results_file));
        rs = results.predictions;
        ds = rs;
end
delete results_file;
switch MODEFL
    case 'classification'

        rs = sign(ds);
    case 'regression'
        rs = ds; 
end

